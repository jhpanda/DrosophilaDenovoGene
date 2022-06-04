""" 
  The script can be used to:
      1. FIND ortholog mapping between species_1 and species_2 from cactus alignment
      2. FIND duplication events from species_1 to species_2 and vice versa

  -input: 
      1. Cactus alignment .hal file, containing genomes of interest
      2. gene based BED files of the two genomes
      3. Aggregated Liftover BED files of the two genomes

  -output: 
      1. Ortholog_mapping
      2. Species-specific duplications
"""

from multiprocessing import Pool
import os,sys,subprocess

class BedLine():
    """ In here, bed lines should contain at least 6 fields
    """
    def __init__(self,line):
        bedinfo = line.strip().split()
        ninfo = len(bedinfo)
        if ninfo>4:
            valid =True
        else:
            valid =False

        try:
            chrom   = bedinfo[0]
        except IndexError:
            chrom   = -1
        try:
            start   = int(bedinfo[1])
        except IndexError:
            start   = -1
        try:
            end     = int(bedinfo[2])
        except IndexError:
            end     = -1
        try:
            name    = bedinfo[3]
        except IndexError:
            name    = "NA"
        try:
            score   = bedinfo[4]
        except IndexError:
            score   = 0
        try:
            strand  = bedinfo[5]
        except IndexError:
            strand  = "."

        others=[] if ninfo<=6 else bedinfo[6:]
        self.chrom = chrom
        self.start = start
        self.end   = end
        self.name  = name
        self.score = score
        self.strand= strand
        self.others= others
        self.valid = valid

def bool_same_access(access1,access2,cutoff):
    """ Make sure that access2 follows access1 in BED file
    """

    flag = 0
    if access1.name == access2.name and access1.chrom == access2.chrom and access1.strand == access1.strand:
        ## we cannot know if the accessions are in ascending or descending orders
        #if abs(access2.start-access1.end)<=cutoff or abs(access1.start-access2.end)<=cutoff:
        #if access2.start-access1.end<=cutoff and access1.start-access2.end<=cutoff:
        bound1 = max(access1.start-cutoff,access2.start-cutoff)
        bound2 = min(access1.end+cutoff,access2.end+cutoff)

        if bound1<bound2:
            flag = 1

    if flag == 1:
        return True
    else:
        return False

def bool_valid_overlap(bedi,bedj,cutoff=0.3):
    """ bedi and bedj are overlapped, this is to determine if they are truly homologs or orthologs, only by percent of matches """
    if bedi.chrom==bedj.chrom and bedi.strand==bedj.strand and bedi.start<bedj.end and bedj.start<bedi.end:
        boundary = [bedi.start,bedi.end,bedj.start,bedj.end]
        boundary.sort() 
        lcommon  = boundary[2] - boundary[1] + 1
        
        #range_i = set(range(bedi.start,bedi.end))
        #range_j = set(range(bedj.start,bedj.end))
 
        #pcommon = range_i.intersection(range_j)
        #lcommon = len(pcommon)
        percent_i = lcommon/(bedi.end-bedi.start+1)
        percent_j = lcommon/(bedj.end-bedj.start+1)
 
        if percent_i>=cutoff and percent_j>=cutoff:
            return True

    return False

def overlap(bedi,bedj):
    lcommon = 0
    if bedi.chrom==bedj.chrom and bedi.strand==bedj.strand and bedi.start<bedj.end and bedj.start<bedi.end:
        boundary = [bedi.start,bedi.end,bedj.start,bedj.end]
        boundary.sort() 
        lcommon  = boundary[2] - boundary[1] + 1
        
    return lcommon

def parseBedLines(fbed,species,abbreviation=True):
    """ A dictionary to store annotation information
    """

    if abbreviation:
        s1i,s1j  = species.split("_")
        s_abbrev = s1i[0].lower()+s1j[0:3]

    bedLines = dict()
    
    namelist = dict()
    lines  = open(fbed,"r")
    for line in lines:
        bedi = BedLine(line)
        name = "%s|%s"%(s_abbrev,bedi.name)
        if name in namelist:
            bedLines[name].append(bedi)
        else:
            bedLines[name]  = [bedi]
            namelist[name] = 1
    lines.close()

    return bedLines

def parseBlastp(fblastp,species1,species2,abbreviation=True):
    """ blastp results are converted to gene pairs
        An example of fblastp contents:
dmag|LOC116936589 dana|LOC6500706
dmag|LOC116936745 dmag|LOC116936745
dmag|YC34_gp10 dpse|HCS66_mgp13;dsim|ATP8;dyak|ATP8;blat|AZ357_gp10;dmel|FBgn0013673;bdor|ATP8;agam|ATP8;amel|ATP8
        gene1 ";"splitted gene2 list
    """
    orthopairs = dict()
    lines = open(fblastp,"r")
    for line in lines:
        if line.startswith(species1) or line.startswith(species2):
            gene1,gene2s = line[:-1].split()
         
            #g1  = genemapping[genes1]
            #g2  = genemapping[genes2]
            orthopairs[gene1] = gene2s.split(";")

    lines.close()
    return orthopairs

def parseSimilarSeq(fsimilar,species1,species2,abbreviation=True):
    """ parse similar sequence to potential ortholog pairs
    """
    orthopairs = dict()
    lines = open(fsimilar,"r")
    for line in lines:
        if line.startswith(species1) or line.startswith(species2):
            g1,g2 = line.split()[0:2]
            ## add g1 and g2 to orthopairs
            if g1 not in orthopairs:
                orthopairs[g1] = [g2]
            else:
                if g2 not in orthopairs[g1]:
                    orthopairs[g1].append(g2)

            if g2 not in orthopairs:
                orthopairs[g2] = [g1]
            else:
                if g1 not in orthopairs[g2]:
                    orthopairs[g2].append(g1)

    lines.close()
    return orthopairs
    

def mapping_one(cds2,genome2,g_liftover1,c_liftover1,species2,orthoMCL_pairs,overlap_cutoff):

    cactus_orthologs  = dict()

    for gene2 in g_liftover1:
        ## only genes that have liftover hits are considered
        #print(gene2,gene2 in orthoMCL_pairs)

        if gene2 not in orthoMCL_pairs:
            for bedi in g_liftover1[gene2]:
                ## if a target gene doesn't have homologs
                ## consider any hit to be unannotated hit, 
                ## even if it is annotated
                ## examples are some overlapping genes
                chrom = bedi.chrom + bedi.strand
                start = str(bedi.start)
                end   = str(bedi.end)
                coor = species2+"|"+":".join([bedi.chrom+bedi.strand,start,end])
                #if species2 == "dyak" and gene2=="dmel|NP_001369097.1":
                #    print(gene2,coor)
                bedi_len = bedi.end-bedi.start+1
                if bedi_len > 0:
                ## of course the length should be larger than 0
                    if gene2 in cactus_orthologs:
                        cactus_orthologs[gene2].append(coor)
                    else:
                        cactus_orthologs[gene2] = [coor]

        elif gene2 in orthoMCL_pairs and gene2 in c_liftover1:
            pairs  = orthoMCL_pairs[gene2]
            pairs  = [g for g in pairs if g.startswith(species2)]
            lmatch = dict()
            lmax   = 0
            for gene1 in pairs:
                ## if a target gene has homologs
                ## in order to be further considered as orthologs
                ## the homolog must matched to it's coding region
                gene1_match_cds = 0
                lcommon = 0
                for bedi in c_liftover1[gene2]:
                    for bedj in cds2[gene1]:
                        lcommon += overlap(bedi,bedj)
                
                if lcommon > lmax:
                    lmax = lcommon
                lmatch[gene1] = lcommon
                if gene2 == "dpse|LOC4803970":
                    print(gene2,gene1,lcommon)

            if lmax:
                for gene1 in lmatch:
                    if lmatch[gene1]/lmax > overlap_cutoff:
                        if gene2 in cactus_orthologs:
                            if gene1 not in cactus_orthologs[gene2]:
                                cactus_orthologs[gene2].append(gene1)
                        else:
                            cactus_orthologs[gene2] = [gene1]
            else:
                ## if the gene does not match to an annotated homolog,
                ## the unannotated hit is not recorded
                ## as there might be duplication, transposition, etc. 
                pass

        elif gene2 not in c_liftover1:
            ## if the hit is not in coding region
            ## might only be overlapping homologs but not paralogs
            ## ignor the hits, neither orthologs nor unannotated hits
            #for bedi in g_liftover1[gene2]:
            #    chrom = bedi.chrom + bedi.strand
            #    start = str(bedi.start)
            #    end   = str(bedi.end)
            #    coor = species2+"|"+":".join([bedi.chrom+bedi.strand,start,end])
            pass

    return cactus_orthologs

#def mapping(species1,species2,dir_orthoMCL,abbreviation=True):
def mapping(species1,species2,fblastp,overlap_cutoff=0.05,abbreviation=True):
    #fbed1 = "../../../2.1_gene_bed_halLiftover/gene_bed/%s-genome.bed"%species1
    #fbed2 = "../../../2.1_gene_bed_halLiftover/gene_bed/%s-genome.bed"%species2
    #fbed3 = "../../../4.1_ortholog_mapping/aggregated_%s_over_%s.bed"%(species2,species1)
    #fbed4 = "../../../4.1_ortholog_mapping/aggregated_%s_over_%s.bed"%(species1,species2)
    fbed1 = "../2.1_gene_bed_halLiftover/gene_bed/%s-genome.bed"%species1
    fbed2 = "../2.1_gene_bed_halLiftover/gene_bed/%s-genome.bed"%species2
    fbed3 = "../1_input_bed_files/%s-genome.bed"%species1
    fbed4 = "../1_input_bed_files/%s-genome.bed"%species2
    fbed5 = "../4.1_ortholog_mapping/aggregated_%s_over_%s.bed"%(species2,species1)
    fbed6 = "../4.1_ortholog_mapping/aggregated_%s_over_%s.bed"%(species1,species2)
    fbed7 = "../4_ortholog_mapping/aggregated_%s_over_%s.bed"%(species2,species1)
    fbed8 = "../4_ortholog_mapping/aggregated_%s_over_%s.bed"%(species1,species2)
    #fbed3 = "../../2.1_gene_bed_halLiftover/%s_over_%s.bed"%(species2,species1)
    #fbed4 = "../../2.1_gene_bed_halLiftover/%s_over_%s.bed"%(species1,species2)

    genome1   = parseBedLines(fbed1,species1,abbreviation)
    genome2   = parseBedLines(fbed2,species2,abbreviation)

    cds1      = parseBedLines(fbed3,species1,abbreviation)
    cds2      = parseBedLines(fbed4,species2,abbreviation)

    g_liftover1 = parseBedLines(fbed5,species1,abbreviation)
    g_liftover2 = parseBedLines(fbed6,species2,abbreviation)

    c_liftover1 = parseBedLines(fbed7,species1,abbreviation)
    c_liftover2 = parseBedLines(fbed8,species2,abbreviation)

    if abbreviation:
        s1i,s1j  = species1.split("_")
        species1 = s1i[0].lower()+s1j[0:3]
        s2i,s2j  = species2.split("_")
        species2 = s2i[0].lower()+s2j[0:3]

    s1 = species1
    s2 = species2

    flocal = "similarPair_%s_%s.txt"%(s1,s2)
    cmd    = "grep -e '%s' %s | grep -e '%s' > %s"%(s1,fblastp,s2,flocal)

    output = subprocess.check_output(cmd,shell=True)
    ## read significant blastp hits
    #fpair = "%s/similarSequence_blastp.txt"%dir_orthoMCL
    #orthoMCL_pairs = parseSimilarSeq_onepair(fpair,species1,species2)
    #fgroup = "%s/genegroups.txt"%dir_orthoMCL
    #orthoMCL_pairs = parseOrthoMCLgroups_onepair(fgroup,species1,species2)
    #orthoMCL_pairs = parseBlastp(flocal,species1,species2)
    orthoMCL_pairs = parseSimilarSeq(flocal,species1,species2)

    os.remove(flocal)

    #print(orthoMCL_pairs)

    ## mapping
    cactus_orthologs_1 = mapping_one(cds2,genome2,g_liftover1,c_liftover1,species2,orthoMCL_pairs,overlap_cutoff)
    n = len(cactus_orthologs_1)
    #for g in cactus_orthologs_1:
    #    print(g)
    #    break
    print("## Forward orthologs %s v.s. %s"%(species1,species2),n)

    cactus_orthologs_2 = mapping_one(cds1,genome1,g_liftover2,c_liftover2,species1,orthoMCL_pairs,overlap_cutoff)
    n = len(cactus_orthologs_2)
    #for g in cactus_orthologs_2:
    #    print(g)
    #    break
    print("## Reverse orthologs %s v.s. %s"%(species1,species2),n)

    cactus_orthologs = {**cactus_orthologs_1,**cactus_orthologs_2}
    return [species1,species2,cactus_orthologs]


class OrthologMapping():
    """ class object to do the mapping
        Input:
          - self.species1
          - self.species2
        Default names of bed files:
          - "../1_input_bed_files/{species}.bed"      annotated BED files
          - "../2_direct_output_bed_files/{species1}_over_{species2}.bed"       halLiftover from species1 to species2
    """
    def __init__(self,genomelist,dir_blastp):
        self.dir_orthoMCL = dir_blastp

        species = []
        lines   = open(genomelist,"r")
        for line in lines:
            species.append(line.strip())
        lines.close()
        self.species = species

        ### ortholog pairs from orthoMCL
        #s1 = "Drosophila_melanogaster"
        #s2 = "Drosophila_simulans"
        #f_ortholog = "%s/pairs/orthologs.txt"%self.dir_orthoMCL
        #pairs = parseOrthoMCLgroups_onepair(f_ortholog,s1,s2)
        #
        #f_coortholog = "%s/pairs/coorthologs.txt"%self.dir_orthoMCL
        #copairs = parseOrthoMCLgroups_onepair(f_coortholog,s1,s2)
        #
        #self.orthoMCL_pairs = {**pairs,**copairs}

        #fpair = "%s/similarSequence_blastp.txt"%dir_orthoMCL
        #self.orthoMCL_pairs = parseSimilarSeq_onepair(fpair,s1,s2)
        #self.orthoMCL_pairs = parseSimilarSeq(fpair)

        ## homolog pairs from blastp results
        #fpair = "%s/similarSequence_blastp2.txt_E0.05_cutoff0.5"%dir_blastp
        fpair = "%s/blastpGenePair_E0.05.txt"%dir_blastp
        self.fpair = fpair
        #self.orthoMCL_pairs = parseBlastp(fpair)
        #print(len(self.orthoMCL_pairs))
    

    def run_mapping(self,ncpu=4,cutoff=300,overlap_cutoff=0.15):
        ## species1 to species2
        ## INPUT:
        ##  Annotation BED: ../1_input_bed_files/{species1}.bed
        ##  halLiftover BED: ../2_direct_output_bed_files/{species1}_over_{species2}.bed
        ## OUTPUT:
        ##  orthlogs1:      find genes in species1 that have synteny maps to genes in species2
        ##  duplications1:  find duplications that are specific to species1


        if ncpu>1:
            pool = Pool(ncpu)
        #s1 = "Drosophila_melanogaster"
        #s2 = "Drosophila_simulans"
        nspecies = len(self.species)
        result_all = []
        #nspecies = 2
        for i in range(nspecies-1):
        #for i in range(3,4):
            for j in range(i+1,nspecies):
        #for i in range(1):
            #for j in range(i+1,nspecies):
        #    for j in range(i+1,2):
            #for j in range(10,11):
        #for i in range(3):
        #    for j in range(i+1,3):
        #for s1 in self.species[0:2]:
        #    for s2 in self.species[0:2]:
        #for s1 in self.species:
        #    for s2 in self.species:
        #        if s1!=s2:
        #        i = 0
        #        j = 3
        #        i = 0
        #        j = 16
        #        j = 14
                s1 = self.species[i]
                s2 = self.species[j]
                #s1 = "dmel"
                #s2 = "dsec"
                #s1 = "Drosophila_melanogaster"
                #s2 = "Drosophila_sechellia"
                if ncpu > 1:
                    #pool.apply_async(mapping,args=(s1,s2,self.dir_orthoMCL),
                    pool.apply_async(mapping,args=(s1,s2,self.fpair,
                                        overlap_cutoff),
                                        callback=result_all.append)
                    print("Starting %s v.s. %s..."%(s1,s2))
                else:
                    #result_all.append(mapping(s1,s2,self.dir_orthoMCL))
                    print("Starting %s v.s. %s..."%(s1,s2))
                    result_all.append(mapping(s1,s2,self.fpair,overlap_cutoff))
                    print("Finish %s v.s. %s"%(s1,s2))
        if ncpu>1:
            pool.close()
            pool.join()
        #cactus_orthologs = mapping(s1,s2,self.orthoMCL_pairs)
        #for s1 in species:
        #    for s2 in species:
        #        if s1!=s2:
        #            mapping(species1,species2,orthoMCL_pairs)

        for species1,species2,cactus_orthologs in result_all:
            print(species1,species2)
            f = open("cactus_orthologs_%s_vs_%s.txt"%(species1,species2),"w")

            all_pairs = dict()
            for g1 in cactus_orthologs:
                for g2 in cactus_orthologs[g1]:
                    if g1 not in all_pairs:
                        all_pairs[g1]  = [g2]
                    else:
                        if g2 not in all_pairs[g1]:
                            all_pairs[g1].append(g2)

            ## aggregate again the hits, since the hits are summarized using proteins, genes with alternative proteins are not aggregated, we need to aggregate these hits again 
            for g1 in all_pairs:
                genes = all_pairs[g1]

                annotated_genes  = []
                unannotated_orfs = []
                for g in genes:
                    if "LOC" not in g and "FBgn" not in g and ":" in g:
                        s,coor = g.split("|")
                        try:
                            chrom,start,end = coor.split(":")
                        except ValueError:
                            print(s,coor)
                            sys.exit()
                        ch      = chrom[-1]
                        chrom   = chrom[0:-1]
                        bedline = "%s %s %s %s %d %s"%(chrom,start,end,g1,0,ch)
                        #if g1 == "dmel|FBgn0265577":
                        #    print(bedline)
                        g_bed   = BedLine(bedline) 
                        #if g1 == "dmel|FBgn0265577":
                        #    print(g_bed.chrom,g_bed.strand)
                        unannotated_orfs.append(g_bed)
                    else:
                        annotated_genes.append(g)

                s = genes[0].split("|")[0]
                aggregated = []
                entry0 = BedLine("")
                #if g1 == "dmel|FBgn0265577":
                #    for orf in unannotated_orfs:
                #        print(orf.chrom,orf.start,orf.end,orf.name,orf.strand)
                unannotated_orfs = sorted(unannotated_orfs,key=lambda x:(x.chrom,x.strand,x.start,x.end))
                #if g1 == "dmel|FBgn0265577":
                #    for orf in unannotated_orfs:
                #        print(orf.chrom,orf.start,orf.end,orf.name,orf.strand)
                for entry in unannotated_orfs:
                    if bool_same_access(entry,entry0,cutoff):
                        entry0.end   = max(entry.end,entry0.end)
                        entry0.start = min(entry.start,entry0.start)
                    else:
                        ### minimum non-genic hit length
                        if entry0.valid and abs(entry0.start-entry0.end)>=overlap_cutoff:
                        #if entry0.valid:
                            coori = "%s|%s%s:%s:%s"%(s,entry0.chrom,
                                                entry0.strand,entry0.start,
                                                entry0.end)
                            aggregated.append(coori)
                        entry0 = entry
                if entry0.valid:
                    coori = "%s|%s%s:%s:%s"%(s,entry0.chrom,
                                        entry0.strand,entry0.start,
                                        entry0.end)
                    ### minimum non-genic hit length
                    if coori not in aggregated and abs(entry0.start-entry0.end)>=overlap_cutoff:
                    #if coori not in aggregated:
                        aggregated.append(coori)
                genes = annotated_genes + aggregated
                pair = "%s %s"%(g1," ".join(genes))
                f.write("%s\n"%pair)
            f.close()

if __name__ == "__main__":
    #s1 = "Drosophila_melanogaster"
    #s2 = "Drosophila_simulans"
    
    genomelist = "fullname_hasHit.txt"
    #dir_orthoMCL = "/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Orthologs/orthoMCL/node062/blastp"
    #blastp_dir = "/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Orthologs/orthoMCL/search_using_blastp/"
    similarSeq_dir = "/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Orthologs/orthoMCL/node062/blastp2/"
    #orthomap = OrthologMapping(genomelist,dir_orthoMCL)
    orthomap = OrthologMapping(genomelist,similarSeq_dir)
    orthomap.run_mapping(ncpu=24,overlap_cutoff=0.05)
    #orthomap.run_mapping(ncpu=1,overlap_cutoff=0.05)
