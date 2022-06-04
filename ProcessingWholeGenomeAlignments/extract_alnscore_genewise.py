

from fasta2seq import *
from Bio.SubsMat import MatrixInfo
from multiprocessing import Pool

import sys
import numpy as np
import scipy.stats
import statsmodels.stats.multitest

def cutoffFromRandomSimulations(X,aln_method,stype):
    """ 
        Scores/Similarity from random genewise/spaln simulations 
        Both mean and standard deviation are fitted by two phase decay
        For fitting procedure, see scorefit.py for detail
    """

    if aln_method=="genewise" and stype=="alnscore":
        ## genewise score ##
        Y01       =   10.999999650
        Plateau1  =    0.010964339
        KFast1    =    0.015192716
        KSlow1    =    0.322977380
        SpanFast1 =    0.148062602
        SpanSlow1 =   10.840972709
        ## std ##
        Y02       =    0.210192576
        Plateau2  =    0.000000102
        KFast2    =    0.020887105
        KSlow2    =    0.001636796
        SpanFast2 =    0.190068618
        SpanSlow2 =    0.020123856

    elif aln_method=="genewise" and stype=="similarity":
        ## similarity ##
        Y01       =    0.696566717
        Plateau1  =    0.005131241
        KFast1    =    0.040691456
        KSlow1    =    0.003112325
        SpanFast1 =    0.642522668
        SpanSlow1 =    0.048912808
        ## std ##
        Y02       =    0.100990127
        Plateau2  =    0.003614450
        KFast2    =    0.012256556
        KSlow2    =    1.718362305
        SpanFast2 =    0.086855728
        SpanSlow2 =    0.010519949

    elif aln_method=="spaln" and stype=="alnscore":
        ## score ##
        Y01       =   17.999999981
        Plateau1  =   -0.039279974
        KFast1    =    0.015644801
        KSlow1    =    0.088245361
        SpanFast1 =    2.736299221
        SpanSlow1 =   15.302980734
        ## std ##
        Y02       =    1.178201033
        Plateau2  =    0.000000000
        KFast2    =    0.001565254
        KSlow2    =    0.036506722
        SpanFast2 =    0.128773380
        SpanSlow2 =    1.049427652

    elif aln_method=="spaln" and stype=="similarity":
        """
        ## similarity ##
        Y01       =    0.372840146
        Plateau1  =    0.204379935
        KFast1    =    0.028078737
        KSlow1    =    1.435058844
        SpanFast1 =    0.156754833
        SpanSlow1 =    0.011705378
        ## std ##
        Y02       =    0.105824272
        Plateau2  =    0.013214079
        KFast2    =    0.002615298
        KSlow2    =  315.961812781
        SpanFast2 =    0.087804776
        SpanSlow2 =    0.004805417
        """
        ## similarity ##
        Y01       =    0.379849567
        Plateau1  =    0.241677256
        KFast1    =    0.020080991
        KSlow1    =    1.721632239
        SpanFast1 =    0.123792600
        SpanSlow1 =    0.014379712
        ## std ##
        Y02       =    2.512711088
        Plateau2  =    0.008828844
        KFast2    =    0.325373367
        KSlow2    =    0.002048836
        SpanFast2 =    2.438368812
        SpanSlow2 =    0.065513433

    expFast1  = np.exp(-KFast1*X)
    expSlow1  = np.exp(-KSlow1*X)
    expFast2  = np.exp(-KFast2*X)
    expSlow2  = np.exp(-KSlow2*X)
    Y         = Plateau1 + SpanFast1*expFast1 + SpanSlow1*expSlow1
    dY        = Plateau2 + SpanFast2*expFast2 + SpanSlow2*expSlow2
    return Y,dY

## probability being de novo of provided similarity and protein length
def probability_denovo(rand_mean,rand_std,score):
    prob = 1-scipy.stats.norm(rand_mean,rand_std).cdf(score)
    return prob

## similar groups according to 
## https://www.bioinformatics.org/sms2/ident_sim.html
similar_grps = "GAVLI FYW CM ST KRH DENQ P X".split()
similar_aa   = dict()
for grp in similar_grps:
    for a in grp:
        similar_aa[a] = grp
#for a in similar_aa:
#    print(a,similar_aa[a])

## use blosum62 as similarity matrix
#  similar_score (or P(i|j)) = 2^((Mij-Mii)/2)
blosum = MatrixInfo.blosum62
aminoacids = "ACDEFGHIKLMNPQRSTVWY"

## lineage groups
def get_lineage(flineage="branch_len.txt"):
    L = "dmel L0 L1 L2 L3 L4 L5 L6 L7 L8 L9".split()
    lines   = open(flineage,"r")
    lineage = dict()
    idx     = 0
    for line in lines:
        l    = L[idx]
        info = line.split()[0].split(",")
        lineage[l] = info
        idx += 1
    return lineage

## start residue
def get_startCodon(seq):
    for s in seq:
        if s in aminoacids:
            return s
    return ""

## add frame shift info to similarity
class AlnScore():
    def __init__(self,score,sim,shift,startCodon,key):
        self.score      = score
        self.shift      = shift
        self.startCodon = startCodon
        self.sim        = sim
        self.key        = key

## similarity between two sequences
def similarity_seq(seq1,seq2,prolen,sim_method="blosum62"):
    slen1 = len(seq1)
    slen2 = len(seq2)
    if slen1!=slen2:
        sim = -1
    else:
        nsim   = 0
        for s1,s2 in zip(seq1,seq2):
            if s1!="-" and s2!="-":
                if sim_method=="group":
                    if s1 in similar_aa[s2]:
                        nsim += 1
                elif sim_method=="blosum62":
                    if s1==s2:
                        nsim += 1
                    else:
                        Mii   = blosum[(s1,s1)]
                        try:
                            Mij   = blosum[(s1,s2)]
                        except KeyError:
                            Mij   = blosum[(s2,s1)]
                        nsim += 2**((Mij-Mii)/2)
                elif sim_method=="ident":
                    if s1==s2:
                        nsim += 1
        sim = min(1,nsim/prolen)
    return sim

def similarity_msa(inp,lineage,aln_method,sim_method,stype,debug):
    """ 
        Infer similarity of unannotated hits by scores or sequence similarity 
        from spliced alignments.
         - aln_method: genewise/spaln
         - sim_method: blosum/group/ident
         - stype:      score tye:   alnscore/similarity
    """
    sequence = fasta2seq(inp)
    keys = [s for s in sequence]
    nkey = len(keys)

    #print(inp)
    if keys==[]:
        print(inp)
    key0 = keys[0]

    ## actual protein length
    #if stype=="alnscore":
    #    prolen = 0
    #    for s in sequence[key0]:
    #        if s!="-":
    #            prolen += 1
    #    prolen = int(prolen/3) -1
    #elif stype=="similarity":
    #    prolen = 0
    #    for s in sequence[key0]:
    #        if s!="-":
    #            prolen += 1
    #    prolen = int(prolen)
    prolen = 0
    for s in sequence[key0]:
        if s!="-":
            prolen += 1
    prolen = int(prolen)

    ## lineage groups in this alignment
    lgrps   = dict()
    species = dict()
    for l in lineage:
        lgrpi = []
        for s in lineage[l]:
            for key in keys:
                if s in key:
                    species[s] = l
                    if s not in lgrpi:
                        lgrpi.append(s)     
        if lgrpi:
            lgrps[l] = lgrpi

    ## unannotated outgroups
    allgrps   = [l for l in lgrps]
    annotated = dict()
    shifted   = dict()
    for l in lgrps:
        annotated[l] = []
        shifted[l]   = []
    for key in sequence:
        s = key[0:4]
        l = species[s]
        if ":" in key:
            annotated[l].append(0)
        else:
            annotated[l].append(1)

        if "shift_Yes" in key:
            shifted[l].append(1)
        else:
            shifted[l].append(0)

    unannotated_grps = []
    shifted_grps     = []
    for l in lgrps:
        if sum(annotated[l]) == 0:
            
            unannotated_grps.append(l)
        if sum(shifted[l]) > 0:
            shifted_grps.append(l)
        
    unannotated_grps = unannotated_grps[::-1]

    ## now pairwise similarity
    if debug:
        print("Detail of pariwsie alnscore between species")
    alnscore = dict()
    for key1 in sequence:
        s1     = key1[0:4]
        seq1   = sequence[key1]
        shift1 = "Yes" if "Yes" in key1 else "No"
        ## start codon of the sequence
        atg1   = get_startCodon(seq1)
        #print(key1)
        #score  = float(key1.split("|")[2][6:])

        shift  = "No" if shift1=="No" else "Yes"
        for key2 in sequence:
            if key2[0:4] == "dmel":

                s2    = key2[0:4]
                pair  = (s1,s2)
                seq2  = sequence[key2]
                atg2  = get_startCodon(seq2)

                score = float(key1.split("|")[2][6:])/prolen
                sim   = similarity_seq(seq1,seq2,prolen,sim_method)

                if debug:
                    print("Protein length %d AlnScore %.3f"%(prolen,score))
                    print(">%s:%s"%(s1,seq1))
                    print(">%s:%s"%(s2,seq2))
                    print("")

                aln1  = AlnScore(score,sim,shift,atg1,key1)
                aln2  = AlnScore(score,sim,shift,atg2,key1)

                if pair in alnscore:
                    if stype == "alnscore":
                        if score>alnscore[pair].score:
                            alnscore[(s1,s2)] = aln1
                            alnscore[(s2,s1)] = aln2
                    elif stype == "similarity":
                        if sim>alnscore[pair].sim:
                            alnscore[(s1,s2)] = aln1
                            alnscore[(s2,s1)] = aln2
                    else:
                        sys.exit("stype can only be alnscore or similarity!")
                else:
                    alnscore[(s1,s2)] = aln1
                    alnscore[(s2,s1)] = aln2


    if debug:
        print("Summary of pariwsie alnscore between species")
        for pair in alnscore:
            print(pair,"%.3f"%alnscore[pair].score)
        print("")

    ## similarity and infer de novo origination
    sum_alnscore   = dict()
    sum_similarity = dict()
    sum_shift      = dict()
    sum_startCodon = dict()
    sum_maxkey     = dict()
    
    info = []
    for l in lgrps:
        maxscore = -1E10
        maxsim   = -1E10
        maxkey   = ""
        maxsps   = ""
        maxlin   = ""
        shift    = ""
        atg      = ""
                            
        for s1 in lgrps[l]:
            s2      = "dmel"
            score_i = alnscore[(s1,s2)].score
            shift_i = alnscore[(s1,s2)].shift
            atg_i   = alnscore[(s1,s2)].startCodon
            sim_i   = alnscore[(s1,s2)].sim
            key_i   = alnscore[(s1,s2)].key
            #print(s1,s2,sim)
            if stype == "alnscore":
                if score_i>maxscore:
                    maxscore = score_i
                    maxkey   = key_i
                    maxsim   = sim_i
                    maxsps   = s2
                    maxlin   = species[s2]
                    shift    = shift_i
                    atg      = atg_i
            if stype == "similarity":
                if sim_i>maxsim:
                    maxscore = score_i
                    maxkey   = key_i
                    maxsim   = sim_i
                    maxsps   = s2
                    maxlin   = species[s2]
                    shift    = shift_i
                    atg      = atg_i
                    score    = score_i

        if l == "dmel" and maxscore == 1E10:
            maxscore = 1E10
            maxsim   = 1.0
            maxlin   = "dmel"
            maxkey   = key0

        #print(l,maxsim,shift,atg)
        sum_alnscore[l]   = maxscore
        sum_similarity[l] = maxsim
        sum_shift[l]      = shift
        sum_startCodon[l] = atg
        sum_maxkey[l]     = maxkey

        annotated = "Unannotated" if l in unannotated_grps else "Annotated"
        if annotated == "Annotated":
            info_s1 = "%s|%s|%s|shift_%s|start_%s|1E10|%.2f"%(l,maxlin,
                                                  annotated,shift,atg,maxsim)
        else:
            info_s1 = "%s|%s|%s|shift_%s|start_%s|%.1f|%.2f"%(l,maxlin,
                                                  annotated,shift,atg,
                                                  maxscore,maxsim)
        info.append(info_s1)

    if debug:
        print("Summary of alnscore in each lineage")
        for l in allgrps:
            print("%-5s %.3f"%(l,sum_alnscore[l]))
        print("")
    # de novo origination
    origin_l = ""
    l0       = ""

    ## mean and std of random simulations ##
    rand_mean,rand_std = cutoffFromRandomSimulations(prolen,aln_method,stype)
    #cutoff = rand_mean + 1.644853626951472*rand_std
    #cutoff = rand_mean + 3.7885560493*rand_std
    #cutoff = rand_mean + 3.09*rand_std

    sigm = 4.75 if fdr_level==1e-6 else 3.09
    cutoff = rand_mean +sigm*rand_std
    #for l in unannotated_grps:
    for l in allgrps[::-1]:
        #if l not in unannotated_grps:
        #    l0 = l
        #    break
        #elif sum_alnscore[l]>0:
        if sum_alnscore[l]>cutoff:
            l0 = l
            break
        else:
            origin_l = l

## print all identified unannotated homolgs
#    for l in allgrps[::-1]:
#        if sum_alnscore[l]>cutoff:
#            if ":" in sum_maxkey[l]:
#                print(key0[0:16],sum_maxkey[l])
## print all identified non-genic
#    for l in allgrps[::-1]:
#        if sum_alnscore[l]<cutoff:
#            if ":" in sum_maxkey[l]:
#                print(key0[0:16],sum_maxkey[l])

    g = key0[0:16]+"|%d"%prolen

    if debug:
        print("Summary:")
        print("Most distant lineage %s; Origination %s"%(l0,origin_l))
        print("%s cutoff %.3f (%.3f+/-%.3f)"%(g,cutoff,rand_mean,rand_std))
        print("")

    #print(g,origin_l,l0)
    if origin_l:
        print(key0[0:16],origin_l,sum_maxkey[origin_l])
        score = sum_alnscore[origin_l]
        sim   = sum_similarity[origin_l]
        prob  = probability_denovo(rand_mean,rand_std,score)
        #print(prolen,score,rand_mean,rand_std)
        info.append("%3.2e"%prob)
        description = ""
        if sum_shift[l0]=="Yes" and sum_startCodon[l0]!="M":
            description = "P_Origin_%s_Shifted&noATG_%s"%(origin_l,l0)
        elif sum_shift[l0]=="Yes" and sum_startCodon[l0]=="M":
            description = "P_Origin_%s_Shifted_%s"%(origin_l,l0)
        elif sum_shift[l0]=="No" and sum_startCodon[l0]!="M":
            description = "P_Origin_%s_noATG_%s"%(origin_l,l0)
        elif sum_shift[l0]=="No" and sum_startCodon[l0]=="M":
            description = "P_Origin_%s"%(origin_l)
        else:
            print(">Warn: no description found for %s!"%g)
    else:
        print(key0[0:16],l0,sum_maxkey[l0])
        score = sum_alnscore[l0]
        sim   = sum_similarity[l0]
        prob  = probability_denovo(rand_mean,rand_std,score)
        #print(prolen,score,rand_mean,rand_std)
        info.append("%3.2e"%prob)
        if sum_shift[l0]=="Yes" and sum_startCodon[l0]!="M":
            description = "P_FrameShift&NoStartCodon_%s"%l0
        elif sum_shift[l0]=="Yes" and sum_startCodon[l0]=="M":
            description = "P_FrameShift_%s"%l0
        elif sum_shift[l0]=="No" and sum_startCodon[l0]!="M":
            description = "P_NoStartCodon_%s"%l0
        elif sum_shift[l0]=="No" and sum_startCodon[l0]=="M":
            description = "P_NotDenovo"
        else:
            print(">Warn: no description found for %s!"%g)
        
    info = " ".join(info)
    result = [g,info,description,prob,sim,score,prolen]
    #print(g,info,description)
    return result


def main():
    inp         = sys.argv[1]
    sim_method  = sys.argv[2]
    aln_method  = sys.argv[3]
    stype       = sys.argv[4]
    sim_methods = "group blosum62 ident".split()
    if sim_method not in sim_methods:
        print("method must be in: ",sim_methods)
        sys.exit(0)
    lineage = get_lineage()
    similarity_msa(inp,lineage,aln_method,sim_method,stype)

def parallel(fgenelist,seq_dir,output_dir="./",aln_method="spaln",sim_method="blosum62",stype="alnscore",fdr_level=1e-3,debug=False,ncpu=30):

    lineage = get_lineage()
    if ncpu>1:
        pool = Pool(ncpu)

    results = []
    lines   = open(fgenelist,"r")
    aln     = "genewise623" if aln_method=="genewise" else "spaln"
    for line in lines:
        pref = line[0:16]
        if "|" in pref:
            pref = pref.replace("|","_")
        inp  ="%s/%s_%s_pro.fa.mafft"%(seq_dir,pref,aln)
        #if stype == "similarity":
        #    inp  ="raw_hit_300_new/%s_%s_pro.fa.mafft"%(pref,aln)
        #else:
        #    inp  ="raw_hit_300_new/%s_%s.fa"%(pref,aln)
        if ncpu>1:
            pool.apply_async(similarity_msa,
                            args=(inp,lineage,aln_method,
                                sim_method,stype,debug),
                            callback=results.append)
        else:
            result = similarity_msa(inp,lineage,aln_method,
                                    sim_method,stype,debug)
            results.append(result)
    if ncpu>1:
        pool.close()
        pool.join()

    prob_list = []
    for result in results:
        g,info,description,prob,sim,score,prolen = result
        prob_list.append(prob)
    reject,Pcorrected = statsmodels.stats.multitest.fdrcorrection(prob_list)

    f1 = open("%s/pro_%s_%s.txt"%(output_dir,stype,aln_method),"w")
    f2 = open("%s/pro_sim_%s_notdenovo.txt"%(output_dir,aln_method),"w")
    f3 = open("%s/pro_sim_%s_denovo.txt"%(output_dir,aln_method),"w")
    f2.write("gene length alnscore similarity prob\n")
    f3.write("gene length alnscore similarity prob\n")
    for result,p in zip(results,Pcorrected):
        g,info,description,prob,sim,score,prolen = result
        if p<fdr_level:
            description_fdr = "Corrected_NotDeNovo"
            f2.write("%s %-5d %5.3f %.3f %3.2e\n"%(g,prolen,score,sim,p))
        else:
            description_fdr = "Corrected_DeNovo"
            f3.write("%s %-5d %5.3f %.3f %3.2e\n"%(g,prolen,score,sim,p))
            
        f1.write("%s %s %s %3.2e %s\n"%(g,info,description,p,description_fdr))
    f1.close()
    f2.close()
    f3.close()


if __name__ == "__main__":
    #main()
    fgenelist = sys.argv[1]
    seq_dir   = sys.argv[2]
    suffix    = sys.argv[3]
    #fgenelist = "annotated_groups.list"
    #fgenelist  = "last_orf_unannotated_nohomologs_v6.txt"
    #fgenelist  = "test.list"
    #aln_method = "genewise"
    #aln_method = "spaln"
    sim_method = "blosum62"
    stype      = "alnscore"
    #stype      = "similarity"
    #fdr_level  = 1e-6
    #output_dir = "./summary_fdr000001/"
    #fdr_levels = [1e3,1e-6]
    #out_dirs   = ["summary_fdr001_%s"%suffix,"summary_fdr000001_%s"%suffix]
    fdr_levels = [1e-6]
    out_dirs   = ["summary_fdr000001_%s"%suffix]
    for fdr_level,output_dir in zip(fdr_levels,out_dirs):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        #debug      = True
        debug      = False
        aln_method = "genewise"
        parallel(fgenelist,seq_dir,output_dir,aln_method,sim_method,stype,fdr_level,debug,ncpu=1)
        aln_method = "spaln"
        parallel(fgenelist,seq_dir,output_dir,aln_method,sim_method,stype,fdr_level,debug,ncpu=1)
