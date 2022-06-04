
lines = open("last_orf_unannotated_v7.txt","r")
ortho_last_orf = dict()
for line in lines:
    info = line.split()
    ortho_last_orf[info[0]] = info
lines.close()

## lineage information ##
lineage = []
lines   = open("branch_len.txt.more","r")
for line in lines:
    grp,t = line[:-1].split()
    grp = grp.split(",")
    lineage += [grp]
lines.close()
lineage = lineage[::-1]

## get the most outgroup species that don't have hits in cactus ##
## we will check if homologs are in these outgroup species ##
## if yes, then this gene is probably not de novo gene ##
cactus_outspecies = dict()
lines = open("last_orf_unannotated_v7.txt","r")
for line in lines:
    fbid       = line[0:16]
    outspecies = []

    nongenic_species = line.split()[-1][0:4]
    for grp in lineage:
        if nongenic_species not in grp:
            for g in grp: 
                outspecies.append(g[0:4])
        else:
            for g in grp:
                if nongenic_species not in g:
                    outspecies.append(g[0:4])
            break
    cactus_outspecies[fbid] = outspecies
lines.close()

## now check if these genes have homolog sequences in outspecies ##
outspecies = "mdom aaeg agam amel aflo mpha bmor dmag".split()
homolog_genes = dict()
lines = open("/ru-auth/local/home/jpeng/scratch/denovogene/revisit/Orthologs/orthoMCL/search_using_blastp/dmel_similarGenes_E0.001.txt","r")
for line in lines:
    fbid       = line[0:16]
    if fbid in ortho_last_orf:
        outspecies = cactus_outspecies[fbid]
        if fbid == "dmel|FBgn0264748":
            print(outspecies)
        if any([s in line for s in outspecies]):
            homolog_genes[fbid] = line[0:-1]
lines.close()

#for key in homolog_genes:
#    print(key,homolog_genes[key])
#key = "dmel|FBgn0262823"
#print(key,homolog_genes[key])

#print(ortho_last_orf["dmel|FBgn00"])
print(ortho_last_orf["dmel|FBgn0264748"])
#for fbid in ["dmel|FBgn0028583"]:
dmel_genes = []
for fbid in ortho_last_orf:
    has_homolog = 0

    for g in ortho_last_orf[fbid]:
        if "dmel" in g and ":" not in g:
            dmel_genes.append(g)

    has_homolog = [g in homolog_genes for g in ortho_last_orf[fbid]]
    has_homolog = any([h for h in has_homolog])
    #print(fbid,has_homolog)
    if not has_homolog:
        #print(" ".join(ortho_last_orf[fbid]))
        pass

#for fbid in dmel_genes:
#    print(fbid)
