"""
    The script can generate ortholog/paralog mapping providing cactus alignments
    Be aware that halLiftover is not capabale of finding duplications of species-specific genes 
"""

import sys

def bool_same_gene(coor_genei,coor_genej,cutoff):
    try:
        si,coori = coor_genei.split("|")
        sj,coorj = coor_genej.split("|")
    except ValueError:
        print(coor_genei,coor_genej)
        sys.exit()

    if coor_genei==coor_genej:
        return coor_genei

    else:
        flag = 0
        if si==sj and ":" in coor_genei and ":" in coor_genej:
            chromi,starti,endi = coori.split(":")
            strandi = chromi[-1]
            chromi  = chromi[0:-1]
 
            chromj,startj,endj = coorj.split(":")
            strandj = chromj[-1]
            chromj  = chromj[0:-1]
 
            if chromi==chromj and strandi==strandj:
 
                starti = int(starti)
                startj = int(startj)
                endi   = int(endi)  
                endj   = int(endj)  
      
                bound1 = max(starti,startj)-cutoff
                bound2 = min(endi,endj)    +cutoff
 
                if bound1 < bound2:
                    return "%s|%s%s:%d:%d"%(si,chromi,strandi,min(starti,startj),max(endi,endj))
    return False

class OrthoMap():
    def __init__(self,gene,cutoff=300):
        self.gene      = gene
        self.paralogs  = []
        self.orthologs = []
        self.has_dup   = False ## if this gene has a duplication?
        self.unique    = True  ## if this gene is species specific?
        self.cutoff    = cutoff

    def add_paralog(self,g):
        self.paralogs.append(g)

    def add_ortholog(self,g):
        self.orthologs.append(g)

def get_species():
    species = dict()
    lines = open("fullname_hasHit.txt","r")
    #lines = open("fullname_test.txt","r")
    for line in lines:
        full = line.strip()
        abbr = full.split("_")
        abbr = abbr[0][0].lower() + abbr[1][0:3]
        species[full] = abbr
    return species

## this function is not updated
def pairwiseOrthologs(self,s1="dmel",s2="dsim"):
    """ orthomaps between two species """
    orthoMap = dict()
    fmap  = "./pairwise_cactus_ortholog/cactus_orthologs_%s_vs_%s.txt"%(s1,s2)
    lines = open(fmap,"r")
    for line in lines:
        if line.startswith(s2):
            linesplit = line[:-1].split()
            g1 = linesplit[0]
            genes_matched = linesplit[1:]
    lines.close()

    return orthoMap

class OrthologsMap():
    def __init__(self,s1="dmel",cutoff=300):
        """ 
            Attributs:
             - self.s1                  focal species
             - self.specific_genes      genes in focal species
             - self.orthomap            ortholog & paralog info of each gene

            note1: In the current state, I only compare dmel to other species

            note2: Ultimately, I want self.orthomap to contain all genes 
                   in all species
        """
        self.s1 = s1
        self.orthomap   = dict()
        self.cutoff     = cutoff

        species = get_species()
        self.species = [species[s] for s in species]

        self.specific_genes = dict()
        lines = open("../../orthoMCL/mapping_genes2protein/%s.genes.unique.txt"%s1,"r")
        for line in lines:
            g = "%s|"%s1 + line.strip()
            self.specific_genes[g] = 1
            self.orthomap[g] = OrthoMap(g)
        lines.close()

        #self.uniqueGene = dict()

    def pairwiseOrthologs(self,s1="dmel",s2="dsim",pair_dir="./pairwise_cactus_ortholog/"):
        """ 1. Identify duplications in s1 relative to s2
            2. Identify one-to-one ortholog between s1 and s2
        """

        genesets = dict()
        ## first assign all ortholog information
        fmap  = "%s/cactus_orthologs_%s_vs_%s.txt"%(pair_dir,s1,s2)
        lines = open(fmap,"r")
        for line in lines:
            linesplit = line[:-1].split()
            g1 = linesplit[0]
            genes_matched = linesplit[1:]

            if g1 not in genesets:
                genesets[g1] = 1 

            if g1 not in self.orthomap:
                self.orthomap[g1] = OrthoMap(g1)
                self.orthomap[g1].unique  = False

            for g2 in genes_matched:
                if g2 not in genesets:
                    genesets[g2] = 1

                ## add g2 to g1's orthologs
                if g2 not in self.orthomap[g1].orthologs:
                    self.orthomap[g1].add_ortholog(g2)

                ## also add g1 to g2's orthologs
                if g2 not in self.orthomap:
                    self.orthomap[g2] = OrthoMap(g2)
                    self.orthomap[g2].unique  = False
                    self.orthomap[g2].add_ortholog(g1)
                        
                else:
                    if g1 not in self.orthomap[g2].orthologs:
                        self.orthomap[g2].add_ortholog(g1)

            ## process duplications
            n_matched = len(genes_matched)
            if n_matched > 1:
                for i in range(n_matched-1):
                    for j in range(i+1,n_matched):
                        g1 = genes_matched[i]
                        g2 = genes_matched[j]

                        self.orthomap[g1].has_dup = True
                        self.orthomap[g2].has_dup = True

                        ## add g2 to g1's paralogs
                        if g2 not in self.orthomap[g1].paralogs:
                            self.orthomap[g1].add_paralog(g2)

                            ## in the meantime, add each of g1's paraplogs to
                            ## g2, and add g2 to each of g1's paralogs's paralog
                            for g3 in self.orthomap[g1].paralogs:
                                if g3 not in self.orthomap[g2].paralogs:
                                    self.orthomap[g2].add_paralog(g3)
                                if g2 not in self.orthomap[g3].paralogs:
                                    self.orthomap[g3].add_paralog(g2)


                        ## add g1 to g2's paralogs
                        if g1 not in self.orthomap[g2].paralogs:
                            self.orthomap[g2].add_paralog(g1)

                            ## in the meantime, add each of g2's paraplogs to
                            ## g1, and add g1 to each of g2's paralogs's paralog
                            for g3 in self.orthomap[g2].paralogs:
                                if g3 not in self.orthomap[g1].paralogs:
                                    self.orthomap[g1].add_paralog(g3)
                                if g1 not in self.orthomap[g3].paralogs:
                                    self.orthomap[g3].add_paralog(g1)
        lines.close()

    def update(self,pair_dir):
        """ update ortholog/paralog relations, all collapse onto focal species
        """
        ns  = len(self.species)
        #ns  = 14
        for i in range(0,ns-1):
            for j in range(i+1,ns):
                s1 = self.species[i]
                s2 = self.species[j]
                #print(s1,s2)
                self.pairwiseOrthologs(s1,s2,pair_dir)

        g = "dsim|Dsimw501_GD16412"
        print(self.orthomap[g].orthologs)
        g = "dmel|FBgn0029590"
        print(self.orthomap[g].orthologs)
        g = "dmel|FBgn0029589"
        print(self.orthomap[g].orthologs)

        #print(self.orthomap["dmel|FBgn0033372"].orthologs)
        groups    = dict()
        for g in self.specific_genes:
            groups[g] = self.orthomap[g].paralogs
            #for g1 in self.orthomap[g].orthologs+self.orthomap[g].paralogs:
            for g1 in self.orthomap[g].paralogs:
                if g1 in self.orthomap:
                    for g2 in self.orthomap[g1].paralogs:
                        if g2 not in self.orthomap[g].paralogs:
                            self.orthomap[g].paralogs.append(g2)
                        if g2 not in groups[g]:
                            groups[g].append(g2)
                    for g2 in self.orthomap[g1].orthologs:
                        if g2 not in self.orthomap[g].orthologs:
                            self.orthomap[g].orthologs.append(g2)

            for g1 in self.orthomap[g].orthologs:
                flag = 0
                for g0 in groups[g]:
                    combined = bool_same_gene(g0,g1,self.cutoff)
                    if combined:
                        groups[g].remove(g0)
                        groups[g].append(combined)
                        flag = 1
                if flag==0:
                    groups[g].append(g1)
                         
                for g2 in self.orthomap[g1].paralogs:
                    if g2 not in groups[g]:
                        flag = 0
                        for g0 in groups[g]:
                            combined = bool_same_gene(g0,g2,self.cutoff)
                            if combined:
                                groups[g].remove(g0)
                                groups[g].append(combined)
                                flag = 1
                        if flag==0:
                            groups[g].append(g2)


        #for g in self.specific_genes:
        #    print(g," ".join(groups[g]))
        #print("dmel|FBgn0000529",groups["dmel|FBgn0000529"])
        #print("dmel|FBgn0051697",groups["dmel|FBgn0051697"])
        #print("dmel|FBgn0051081",groups["dmel|FBgn0051081"])
        #print("dmel|FBgn0051080",groups["dmel|FBgn0051080"])
        print("dmel|FBgn0029590",groups["dmel|FBgn0029590"])
        print("dmel|FBgn0029589",groups["dmel|FBgn0029589"])
        #g = "dmel|FBgn0051081"
        #print(self.orthomap[g].paralogs)

        ncopy = dict()
        groups_rm = []     ## duplications are included in only one group
        genes_rm  = dict() ## add a duplicate gene to genes_rm, and it will not
                           ## be considered again, this is to get groups_rm
        for s in self.species:
            ncopy[s] = 0
        idx = 0
        for g in self.specific_genes:
            if g not in genes_rm:
                if not self.orthomap[g].has_dup:
                    ncopy[self.s1] = 1
                for g1 in groups[g]:
                    s = g1[0:4]
                    ncopy[s] += 1

                var = [ncopy[s] for s in self.species]
                #print(g,var)

                #if var[1]>var[0]:
                #    print("")
                #    print(g," ".join(groups[g]))
                #    print("")

                #if len(set(var)) != 1:
                #if len(set(var)) != 1 and max(var)>1:
                #if len(set(var)) >= 1:
                    #groups_rm += [groups[g]]
                if self.orthomap[g].has_dup:
                    #print(" ".join(groups_rm[idx]))
                    groups_rm += [groups[g]]
                else:
                    groups_rm += [[g]+groups[g]]
                    #print(g," ".join(groups_rm[idx]))
                idx += 1
                    
                for g1 in self.orthomap[g].paralogs:
                    genes_rm[g1] = 1

                for s in self.species:
                    ncopy[s] = 0

        for grp in groups_rm:
            newgrp = [s for s in grp]
            n_grp  = len(grp)
            for i in range(n_grp-1):
                for j in range(i+1,n_grp):
                    g1 = grp[i]
                    g2 = grp[j]
                    combined = bool_same_gene(g1,g2,self.cutoff)
                    if combined:
                        newgrp[i] = combined
                        newgrp[j] = ""
            newgrp = [s for s in newgrp if s!=""]
            
            for g in newgrp:
                s = g[0:4]
                ncopy[s] += 1
            
            var = [ncopy[s] for s in self.species]
            #if len(set(var)) > 1 and max(var)>1:
            print(" ".join(newgrp))
                    
            for s in self.species:
                ncopy[s] = 0



        #f = open("dmel_orthoMap.txt","w")
        #f.write("Gene\tParalogs\tOrthologs\n")




        #print(self.orthoMap["dmel|FBgn0287183"].paralogs)
        #f.close()

        #f = open("dmel_dup.txt","w")
        #existed = dict()
        ##for g in paralog_groups:
        #    group = paralog_groups[g]
        #    tup   = tuple(group)
        #    if tup not in existed:
        #        f.write("%s\n"%" ".join(group))
        #        existed[tup] = 1
        #f.close()
#
#
if __name__ == "__main__":
    ortho = OrthologsMap()

    pair_dir = "./pairwise_cactus_ortholog/"
    ortho.update(pair_dir) 
