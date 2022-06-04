
import sys

def read_branches():
    f_lineage = "branch_len.txt"
    lineage = []
    distance= []
    lines   = open(f_lineage,"r")
    for line in lines:
        grp,t = line[:-1].split()
        grp = grp.split(",")
        lineage += [grp]
        distance.append(float(t))

    nl = len(lineage)
    print(lineage)
    return lineage


def main(inp,out):

    lineage_groups = read_branches()

    lines = open(inp,"r")
    idx = 0

    f = open(out,"w")
    for line in lines:
        grp = line.strip().split()

        dmel_genes  = []
        other_genes = []
        for g in grp:
            if "dmel" in g:
                dmel_genes.append(g)
            else:
                other_genes.append(g)

        branch_groups = []
        for lgrp in lineage_groups:
            gene_group = []
            for g in grp:
                s = g[0:4]
                if s in lgrp:
                    gene_group.append(g)
            for g in gene_group:
                grp.remove(g) 

            if gene_group:
                branch_groups += [gene_group]

        if branch_groups:
            flag = 1
            for g in branch_groups[-1]:
                if ":" not in g: 
            ## not consider this gene if any gene in last group is annotated
                    flag = 0

            if flag == 1:
                f.write(line)
            #print(line.strip())
            #print(idx,len(grp))
        #idx += 1

    lines.close()
    f.close()

if __name__ == "__main__":
    
    inp = sys.argv[1]
    #inp = "unique_groups.txt"
    out = sys.argv[2]
    #out = "last_orf_uannotated.txt"
    main(inp,out)
