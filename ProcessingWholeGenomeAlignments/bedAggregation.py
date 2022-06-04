
"""
    The script can aggregate many tiny hits from halLiftover.
"""

import os,sys,gzip
from multiprocessing import Pool

class BedLine():
    def __init__(self,line = ""):
        bedinfo = line[0:-1].split()
        ninfo   = len(bedinfo)

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
            strand  = "+"

        if ninfo<4:
            valid = False
        else:
            valid = True
          
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
        ##if abs(access2.start-access1.end)<=cutoff or abs(access1.start-access2.end)<=cutoff:
        if access2.start-access1.end<=cutoff and access1.start-access2.end<=cutoff:
            flag = 1

    if flag == 1:
        return True
    else:
        return False

def bedAggregation(fbed,cutoff):
    name = os.path.basename(fbed)
    if fbed.endswith(".gz"):
        lines = gzip.open(fbed,"rt")
        name = name[0:-7]
    elif fbed.endswith(".bed"):
        lines = open(fbed,"r")
        name = name[0:-4]
    else:
        print("Unknown foramt %s"%fbed)
        sys.exit(0)


    ## first read all entries
    allEntries = dict()
    for line in lines:
        entry_i = BedLine(line)
        name_i  = entry_i.name

        if name_i in allEntries:
            allEntries[name_i].append(entry_i)
        else:
            allEntries[name_i] = [entry_i]

    ## then analyze entries with the same entry name to see if they can be combined 
    f = open("aggregated_%s.bed"%name,"w")
    for name in allEntries:
        if len(allEntries[name]) > 1:
            entries = sorted(allEntries[name],key=lambda x:(x.chrom,x.strand,x.start,x.end))
        else:
            entries = allEntries[name]

        aggregated = []
        entry0 = allEntries[name][0] 
        for entry in entries:
            if bool_same_access(entry,entry0,cutoff):
                entry0.end   = max(entry.end,entry0.end)
                entry0.start = min(entry.start,entry0.start)
            else:
                aggregated.append(entry0)
                #if entry0.valid:
                #    bedinfo  = [entry0.chrom,str(entry0.start),
                #                str(entry0.end),entry0.name,
                #                access_ctl.score,entry0.strand
                #               ]
                #    bedinfo += entry0.others
                #    f.write("%s\n"%("\t".join(bedinfo)))
                entry0 = entry

        flag = 1
        for agg in aggregated:
            if bool_same_access(entry0,agg,cutoff):
                aggregated.remove(agg)
                aggregated.append(entry0)
                flag = 0
        if flag:
            aggregated.append(entry0)

        aggregated = sorted(aggregated,key=lambda x:(x.chrom,x.strand,x.start,x.end))
        existed_bed = []
        for entry0 in aggregated:
            bedinfo  = [entry0.chrom,str(entry0.start),
                        str(entry0.end),entry0.name,
                        entry0.score,entry0.strand
                       ]
            bedinfo += entry0.others
            if bedinfo not in existed_bed:
                #if entry0.end - entry0.start>=30:
                f.write("%s\n"%("\t".join(bedinfo)))
                existed_bed.append(bedinfo)
            #f.write("%s\n"%("\t".join(bedinfo)))
    print("Finished %s"%fbed)
                
    f.close()

def main(flist,ncpu=1,cutoff=100):
    species = []
    lines = open(flist,"r")
    for line in lines:
        species.append(line.strip()) 

    if ncpu>1:
        pool = Pool(ncpu)
    for s1 in species:
        for s2 in species:
            if s1!=s2:
                bed = "../2.1_gene_bed_halLiftover/%s_over_%s.bed"%(s1,s2)
                if ncpu>1:
                    pool.apply_async(bedAggregation,args=(bed,cutoff))
                else:
                    bedAggregation(bed,cutoff)
    if ncpu>1:
        pool.close()
        pool.join()

    
if __name__ == "__main__":

    cutoff = 100

    """
    s1 = "Drosophila_simulans"
    s2 = "Drosophila_melanogaster"

    bedpath = "../2_direct_output_bed_files/"
    bed = "%s/%s_over_%s.bed"%(bedpath,s1,s2)

    agg = BedAggregation()
    agg.aggregate(bed)
    """

    flist  = "fullname_hasHit.txt"
    main(flist,ncpu=20,cutoff=cutoff)
