
import os,sys,subprocess
from multiprocessing import Pool
from translate_cds import *
from fasta2seq import *

### global variables ###
base1 = "A C G T N - R Y S W K M H D B V".split()
base2 = "T G C A N - Y R W S M K D H V B".split()
"""
Some additional bases
R   A or G
Y   C or T
S   G or C
W   A or T
K   G or T
M   A or C
B   C or G or T
D   A or G or T
H   A or C or T
V   A or C or G
"""
basepair = dict(zip(base1,base2))

#extension       = "coverage50_orthoMCL"
#extension       = "blastp_E1.00"
#extension       = "new_v7"
extension       = ""

genome_dir      = "/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Genomes.masked_tmpfiles/"
geneBed_dir     = "/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Orthologs/cactus/2.1_gene_bed_halLiftover/gene_bed/"
proBed_dir      = "/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Orthologs/cactus/1_input_bed_files/originBed/"
proSeq_dir      = "/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Orthologs/cactus/5.2_gene_bed_duplications_match0.5_E0.05/"
geneMapping_dir = "/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Orthologs/orthoMCL/mapping_genes2protein" 
fspecies        = "abbrev.txt"

if extension:
    outdir          = "raw_hit_300_%s"%extension
else:
    outdir          = "raw_hit"

if not os.path.isdir(outdir):
    os.mkdir(outdir)

genewise = "/ru-auth/local/home/lzhao/.linuxbrew/bin/genewise"
aln = "/ru-auth/local/home/jpeng/scratch/softwares/prrnaln5.2.0.linux64/bin/aln"
#spaln = "/ru-auth/local/home/jpeng/scratch/softwares/spaln2.4.4.linux64/bin/spaln"
spaln = "spaln"

def exon_num_from_bed(fbed,species):
    ## get exon numbers of proteins from bed files
    exon_num = dict()
    lines    = open(fbed,"r")
    for line in lines:
        chrom,start,end,pro,score,strand,*_ = line.split()
        pro = "%s|%s"%(species,pro)
        if pro in exon_num:
            exon_num[pro] += 1
        else:
            exon_num[pro]  = 1
    return exon_num

def genewise_comp(ref,key,cds,fout,mode,method):
    if method == "genewise623":
        #sys.stdout.write("processing %s-%s...\n"%(ref,cds))
        #cmd = "%s %s %s -cdna -quiet | awk 'NR==1{printf $0\"\\n\";next}{printf /^>/ ? \"\\n\"$0\"\\t\" : $0}' $tmp  | sed -e 's/\/\///g'  -e 's/-mRNA-1..*sp/-mRNA-1/g' | grep \">\" | awk '{print$1\"\\n\"$2}' | sed -e \"/^>/d\" | awk '{printf(\"%%s\",$0)}'"%(genewise,ref,cds)
        cmd = "%s %s %s -sum -gff"%(genewise,ref,cds)
        print(cmd)
        ## read dna
        dnaseq = fasta2seq(cds)
        for key in dnaseq:
            dnaseq = dnaseq[key]

        output = subprocess.check_output(cmd,shell=True).decode(sys.stdout.encoding)
        output = output.split("\n")
        nlines = len(output)

        outputseq = ""
        nexon     = 0
        score     = 0
        idx       = 0
        shift     = "No"
        idx       = 0
        while idx < nlines:
            if output[idx].startswith("Bits"):
                break
            idx += 1
        idx += 1
        score,query,qstart,qend,target,sstart,send,n_indels,n_introns = output[idx].split()
        score  = float(score)
        nindel = int(n_indels)
        nexon  = int(n_introns) + 1
        
        shift  = "No" if nindel==0 else "Yes"
        #print(score,shift,nexon)

        while idx < nlines:
            if 'cds' in output[idx]:
                chrom,src,stype,start,end,sc,strand,*_ = output[idx].split()
                seq = dnaseq[int(start)-1:int(end)]
                if strand == "-":
                    seq = [basepair[s] for s in seq[::-1]]
                    seq = "".join(seq)
                    outputseq = "%s%s"%(seq,outputseq)
                else:
                    outputseq = "%s%s"%(outputseq,seq)
            idx += 1

        seqline = ">%s|score_%.1f|shift_%s|%d\n%s\n"%(key,score,shift,nexon,outputseq)
        return seqline

    elif method == "spaln":
        ## produce gff file and determine exon numbers
        cmd    = "%s -Q0 -O0 -TInsectDm -yL 45 -pw %s %s"%(spaln,cds,ref)
        output = subprocess.check_output(cmd,stderr=subprocess.STDOUT,
                    shell=True).decode(sys.stdout.encoding)
        print(cmd)

        ## read dna
        dnaseq = fasta2seq(cds)
        for key in dnaseq:
            dnaseq = dnaseq[key]

        ## process output
        #print(output)
        output    = output.split("\n")
        nlines    = len(output)
        outputseq = ""
        nexon     = 0
        score     = 0
        idx       = 0
        shift     = "No"
        while idx < nlines:
            if output[idx] == "##":
                pass
            elif 'Frame' in output[idx]:
                shift = "Yes"
            elif 'gene' in output[idx]:
                chrom,src,stype,start,end,sc,strand,*_ = output[idx].split()
                score = float(sc)
            elif 'cds' in output[idx]:
                chrom,src,stype,start,end,sc,strand,*_ = output[idx].split()
                seq = dnaseq[int(start)-1:int(end)]
                if strand == "-":
                    seq = [basepair[s] for s in seq[::-1]]
                    seq = "".join(seq)
                    outputseq = "%s%s"%(seq,outputseq)
                else:
                    outputseq = "%s%s"%(outputseq,seq)
                nexon += 1
            idx += 1
        seqline = ">%s|score_%.1f|shift_%s|%d\n%s\n"%(key,score,shift,nexon,outputseq)
        #print(seqline)
        return seqline

class LineageSequenceAnalyzer():
    """ Compare genes that have unannotated hit in the most outgroup sequence 
    """
    def __init__(self,outdir,cutoff=150):
        ## a cutoff is set so we can add +/- cutoff base pairs to genes ##
        self.cutoff = cutoff
        self.outdir = outdir
    
        ## read in species to analyze ##
        lines = open(fspecies,"r")
        self.species = []
        for line in lines:
            self.species.append(line.strip())
        lines.close()
        self.nspecies = len(self.species)

        ## read in protein names of genes, further act as genewise reference ##
        self.genemapping = dict()
        for s in self.species:
            fmap = "%s/%s.mapping"%(geneMapping_dir,s)
            lines = open(fmap,"r")
            for line in lines:
                pro,gene = line.strip().split()
                pro  = s+"|"+pro
                gene = s+"|"+gene
                if gene not in self.genemapping:
                    self.genemapping[gene] = [pro]
                else:
                    if pro not in self.genemapping[gene]:
                        self.genemapping[gene].append(pro)
            lines.close()

        ## read in exon numbers for annotated protein-coding genes
        pro_exon_num = dict()
        for s in self.species:
            fbed = "%s/%s-genome.bed"%(proBed_dir,s)
        
            if not os.path.isfile(fbed):
                print("%s not exist!"%fbed)
                sys.exit(0)
        
            lines = open(fbed,"r")
            for line in lines:
                chrom,start,end,pro,score,strand = line[:-1].split()
                pro = "%s|%s"%(s,pro)
                if pro in pro_exon_num:
                    pro_exon_num[pro] += 1
                else:
                    pro_exon_num[pro]  = 1
            lines.close()

        ## assigning pro exon numbers to genes
        gene_exon_num = dict()
        for gene in self.genemapping:
            nexon = 0
            for pro in self.genemapping[gene]:
                if pro in pro_exon_num:
                    if pro_exon_num[pro] > nexon:
                        nexon = pro_exon_num[pro]
            gene_exon_num[gene] = nexon
        self.gene_exon_num = gene_exon_num

        ## read in protein sequences, act as genewise reference ##
        self.cdsSeq = dict()
        for s in self.species:
            fcds = "%s/cds_genomes/%s_cds_sequences.fasta"%(proSeq_dir,s)
            self.cdsSeq = {** self.cdsSeq, ** fasta2seq(fcds)}

        self.proSeq = fasta2seq("%s/pro_genomes/dmel_pro_sequences.fasta"%proSeq_dir)

        print("LineageSequenceAnalyzer initialized!")

    def load_genomes(self):
        ## read in genome sequences ##
        self.genomes = dict()
        for s in self.species:
            fgenome  = "%s/%s-genome.fasta"%(genome_dir,s)
            if os.path.isfile(fgenome):
                self.genomes[s] = fasta2seq(fgenome)
            else:
                print("%s not exist!"%fgenome)
                sys.exit(0)

        ## read in gene bed and extract coordinates of annotated gene ##
        self.geneCoor = dict()
        self.geneLen  = dict()
        for s in self.species:
            fbed = "%s/%s-genome.bed"%(geneBed_dir,s)

            if not os.path.isfile(fbed):
                print("%s not exist!"%fbed)
                sys.exit(0)

            lines = open(fbed,"r")
            for line in lines:
                chrom,start,end,gene,score,strand = line[:-1].split()
                coor = "%s%s:%s:%s"%(chrom,strand,start,end)
                gene = s+"|"+gene
                self.geneCoor[gene] = s+"|"+coor
                self.geneLen[gene]  = int(end)-int(start)+1
            lines.close()

    ## for annotated genes, extract coding sequence
    def extractSeqFromGene(self,gene):
        gene = gene.replace("|","_")
        seq  = self.cdsSeq[gene]
        return seq

    ## for unannotated genes, extract the hits
    def extractSeqFromCoor(self,g_coor,reflen,trb):
        s,coor   = g_coor.split("|")
        sequence = self.genomes[s]
        chrom,start,end = coor.split(":")
        strand = chrom[-1]
        chrom  = chrom[:-1]
        start  = int(start)
        end    = int(end)

        #hit_len = end-start+1
        #trb     = (reflen-hit_len+trb)//2
        #trb     = 0 if trb<0 else trb
 
        start = 0 if start-1-trb<0 else start-1-trb
        end   = end + trb

        try:
            seq = sequence[chrom][start:end]
        except KeyError:
            seq = ""
        if strand == "-":
            seq = [basepair[s] for s in seq[::-1]]
            seq = "".join(seq)
        return seq

    def extractSeqFromGroup(self,group,name):
        group_seq = dict()
        #g0        = group[0]
        #group_seq[g0] = self.cdsSeq[name]
        #for g in group[1:]:
        key    = name.replace("_","|")
        reflen = self.geneLen[key]

        for g in group:
            if ":" in g:
                trb  = self.cutoff
                seq  = self.extractSeqFromCoor(g,reflen,trb)
                key  = g
            else:
                #g_coor = self.geneCoor[g]
                #trb  = 0
                #seq  = self.extractSeqFromCoor(g_coor,trb)
                seq  = self.extractSeqFromGene(g)
                key  = "%s|score_1E10|shift_No|%d"%(g,self.gene_exon_num[g])
            group_seq[key] = seq
        seq2fasta(group_seq,"%s/%s.fasta"%(self.outdir,name))

    def extractGroupList(self,fgroup,genelist):
        lines = open(fgroup,"r")
        idx = 0
        for line in lines:
            fbid    = line[0:16]
            if fbid in genelist:
                group   = line[:-1].split()
                #grpname = "group_%d"%idx
                grpname = fbid.replace("|","_")
                self.extractSeqFromGroup(group,grpname)
                idx += 1

    def genewiseGroup(self,grpname,refgene,method="genewise"):

        #try:
        #    proteins = self.genemapping[refgene]
        #except KeyError:
        #    proteins = []
        #
        #print(proteins)

        

        f2align  = "%s/%s.fasta"%(self.outdir,grpname)
        dnaseq   = fasta2seq(f2align) 

        #proteins = proteins[0:1]
        ndna = len(dnaseq)

        #N = npro*ndna
        #if N:
        #    pool = Pool(N)

        name = refgene.replace("|","_")

        #p_idx = 0
        #for pro in proteins:
        #    proseq = dict()
        ##    proseq[pro] = self.proSeq[pro]
        #    if refgene == "dmel|FBgn0051439":
        #        proseq[pro] = proseq[pro].replace("TTT","T")
        proseq = {name:self.proSeq[name]}

        ref = "%s/%s_temp_pro.fa"%(self.outdir,name)
        seq2fasta(proseq,ref)

        d_idx  = 0
        fout   = "%s/%s_%s.fa"%(self.outdir,grpname,method)

        #ndna = 1
        if ndna>1:
            pool   = Pool(ndna)
        result = []
        for key in dnaseq:
            cds = "%s/%s_temp_dna_%d.fa"%(self.outdir,name,d_idx)
            f = open(cds,"w")
            f.write(">%s\n%s\n"%(key,dnaseq[key]))
            f.close()
        
            if ":" in key:
                mode = "w" if d_idx==0 else "a"
                #pool.apply_async(genewise_comp,args=(ref,key,cds,fout,mode,method))
                if ndna>1:
                    pool.apply_async(genewise_comp,
                                     args=(ref,key,cds,fout,mode,method),
                                     callback=result.append)
                else:
                    r = genewise_comp(ref,key,cds,fout,mode,method)
                    result.append(r)

            else:
                r = ">%s\n%s\n"%(key,dnaseq[key])
                result.append(r)

            d_idx += 1

        if ndna>1:
            pool.close()
            pool.join()
        #print(result)

        result_sort = []

        for r in result:
            if refgene in r:
                result_sort.append(r)
                result.remove(r)

        i = 0
        while i<self.nspecies:
            flag = 0
            for r in result:
                if self.species[i] in r:
                    result_sort.append(r)
                    result.remove(r)
                    flag = 1
            if flag==0:
                i += 1

        with open(fout,"w") as f:
            for seqline in result_sort:
                f.write(seqline)

        #if N:
        #    pool.close()
        #    pool.join()

    def genewiseAll(self,fgroup,genelist,method="genewise623",remove=True):
        lines = open(fgroup,"r")
        idx = 0
        #pool = Pool(ncpu)
        for line in lines:
            fbid    = line[0:16]
            if fbid in genelist:
                group   = line.split()
                grpname = fbid.replace("|","_")
                refgene = group[0]
                print(grpname,refgene)
                self.genewiseGroup(grpname,refgene,method)
                idx += 1

        if remove:
            for f in os.listdir("./"):
                #if "_temp_dna_" in f or "_temp_pro_" in f:
                if "_temp_dna_" in f:
                    os.remove(f)


if __name__ == "__main__":
    #fgroup = "last_orf_unannotated_nohomologs_v7.txt"
    #outdir = "raw_hit_300_new_v7"
    fgroup = sys.argv[1]
    outdir = sys.argv[2]
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    lines = open(fgroup,"r")
    #lines = open("test.list","r")
    genelist = dict()
    for line in lines:
        genelist[line[0:16]] = 1
    #genelist["dmel|FBgn0085361"] = 1
    
    #novoseq = LineageSequenceAnalyzer(outdir,cutoff=300)
    novoseq = LineageSequenceAnalyzer(outdir,cutoff=0)
    novoseq.load_genomes()
    novoseq.extractGroupList(fgroup,genelist)
    #novoseq.genewiseAll(fgroup,genelist,method="genewise623",remove=False)
    #novoseq.genewiseAll(fgroup,genelist,method="spaln",remove=False)
