"""
    This script, I will follow deepMSA.
"""

import os,sys,shutil,subprocess,datetime
import glob
from fasta2seq import *

#evalues=[1e-10,1e-8,1e-6,1e-4,1e-2]
#evalues=[1]
ecut = 0.01

database = {
            "uniclust"  :"/ru-auth/local/home/jpeng/scratch/softwares/deepMSA_hhsuite3/hhsuite2/database/UniRef30_2020_06",
            "uniref"    :"/ru-auth/local/home/jpeng/scratch/softwares/deepMSA_hhsuite3/hhsuite2/database/uniref90.fasta",
            "metaclust" :"/ru-auth/local/home/jpeng/scratch/softwares/deepMSA_hhsuite3/hhsuite2/database/metaclust_nr.fasta",
            "tara"      :"/ru-auth/local/home/jpeng/scratch/softwares/deepMSA_hhsuite3/hhsuite2/database/tara.fasta"
           }

def is_empty_file(inp):
    """ file is empty if its size is 0"""
    return os.path.exists(inp) and os.stat(inp).st_size == 0

def chop_seq(inp):
    sequence = fasta2seq(inp)
    for key in sequence:
        seq_ref = sequence[key]
        break
    seq_len = len(seq_ref)
    sequence_new = dict()
    for key in sequence:
        seq = []
        for i in range(seq_len):
            if seq_ref[i]!='-':
                try:
                    seq.append(sequence[key][i])
                except IndexError:
                    seq.append('-')
        sequence_new[key] = ''.join(seq)
    seq2fasta(sequence_new,inp)

def filterNonRegularResidues(inp):
    sequence = fasta2seq(inp)
    for key in sequence:
        seq = sequence[key]
        seq = seq.replace('U','C')
        seq = seq.replace('u','C')
        sequence[key] = seq
    seq2fasta(sequence,inp)
    

def check_Nf(inp):

    prefix = os.path.splitext(inp)[0]
    a3mout = "%s.checkNf.a3m"%(prefix)
    #print(prefix,a3mout)
    lines = open(inp,'r')

    with open(a3mout,'w') as f:
        for line in lines:
            if not line.startswith('>'):
                f.write(line)
    cmd = "calNf %s"%a3mout
    nf  = subprocess.check_output(cmd,shell=True).decode('utf-8')
    nf  = float(nf)
    return nf
    

def enough_msa(prefix,db):
    inp    = "%s.%s.a3m"%(prefix,db)

    nf = check_Nf(inp)
    time = datetime.datetime.now().time()
    sys.stdout.write("-%s Nf of %s is: %f\n"%(time,inp,nf))
    return True if nf>=64 else False

def mafft_realign(inp,out,ncpu):
    """ Realigning the sequences are critical!!!!!! """

    filterNonRegularResidues(inp)
    seq = fasta2seq(inp)
    cmd = 'mafft-linsi --thread %d %s > %s'%(ncpu,inp,out)
    if len(seq)>200:
        cmd = 'mafft --auto --thread %d %s > %s'%(ncpu,inp,out)
    time = datetime.datetime.now().time()
    sys.stdout.write("-%s ReAlign %s to %s: %s\n"%(time,inp,out,cmd))
    subprocess.check_output(cmd,shell=True)
    chop_seq(out)

def hhblits(inp,a3mout,hhmout,ncpu,db):
    a3m   = a3mout
    hhm   = hhmout
    query = inp
    cmd   = "hhblits -cpu %d -id 99 -diff inf -cov 50 -n 3 -maxmem 16.0 "\
            "-i %s -d %s -maxres 65535 "\
            "-oa3m %s -ohhm %s"%(ncpu,query,database[db],a3m,hhm)

    ## Let's do it using sbatch ##
    sbatch = "sbatch -N 1 -n %d -o hhblits_%s.out --wait --wrap="%(ncpu,db)
    sbatch_cmd = "%s\"%s\""%(sbatch,cmd)

    #print(sbatch_cmd)
    time = datetime.datetime.now().time()
    sys.stdout.write("\n-%s Start HHBLits: %s\n\n"%(time,sbatch_cmd))
    subprocess.check_output(sbatch_cmd,shell=True)
    #if not db == 'uniclust':
    ##subprocess.check_output(cmd,shell=True)

    # realign output MSA ##
    prefix    = os.path.splitext(a3m)[0] 
    a3m_realn = '%s.realn.a3m'%prefix
    a3m_bak   = '%s.bak'%a3m

    mafft_realign(a3m,a3m_realn,ncpu)
    shutil.copy(a3m,a3m_bak)
    shutil.copy(a3m_realn,a3m)

def hhblits_on_db(prefix,ncpu,db,skip_search=False):
    """ db could be custom db built from hmmer search; or uniclust30 """
    inp    = "%s.fasta"%prefix
    a3mout = "%s.%s.a3m"%(prefix,db)
    hhmout = "%s.%s.hhm"%(prefix,db)

    time = datetime.datetime.now().time()
    sys.stdout.write("\n-%s Starting HHBLits on %s \n\n"%(time,db))
    if not skip_search:
        hhblits(inp,a3mout,hhmout,ncpu,db)
    else:
        time = datetime.datetime.now().time()
        sys.stdout.write("\n-%s Skipping HHBLits on %s \n\n"%(time,db))

    ## finally to make sure all alignments have the same colum numbers ##
    #chop_seq(a3mout) ## sometimes hhblits output can be out of format 

def jackhmmer_search(prefix,previous_db,db,skip_search=False):
    query = "%s.fasta"%prefix

    ## search ##
    out = '%s.%s.sto'%(prefix,db)
    tmp = "%s.%s.out"%(prefix,db)
    tbl = "%s.%s.tbl"%(prefix,db)

    cmd = "qjackhmmer --cpu %d -N 3 -A %s -o %s -E 10 --incE 1e-3 "\
          "--tblout %s %s %s"%(ncpu,out,tmp,tbl,query,database[db])

    time = datetime.datetime.now().time()

    ## sbatch 
    sbatch = "sbatch -N 1 -n %d -o jackhmmer_%s.out --wait --wrap="%(ncpu,db)
    sbatch_cmd = "%s\"%s\""%(sbatch,cmd)

    sys.stdout.write("\n-%s qjackhmmer on %s: %s\n\n"%(time,db,cmd))
    if not skip_search:
        #subprocess.check_output(cmd,shell=True)
        subprocess.check_output(sbatch_cmd,shell=True)
    else:
        time = datetime.datetime.now().time()
        sys.stdout.write("\n-%s Skipping qhmmsearch on %s\n\n"%(time,db))

    ## specify output file names ##
    a3mout = '%s.%s.a3m'%(prefix,db)

    ## if no alignments returned ##
    if is_empty_file(out):
        shutil.copy(a3m,a3mout)
        return 0
    else:
        return 1

def hmmer_search(prefix,previous_db,db,skip_search=False):
    query = "%s.fasta"%prefix
    a3m   = "%s.%s.a3m"%(prefix,previous_db)
    hmm   = "%s.%s.hmm"%(prefix,previous_db)
    #chop_seq(a3m) ## sometimes hhblits output can be out of format 
    cmd = "hmmbuild %s %s"%(hmm,a3m)

    time = datetime.datetime.now().time()
    print('\n',"-%s Building HMM profile from previous MSA: "%time,cmd,'\n')
    subprocess.check_output(cmd,shell=True)

    ## then search ##
    out = '%s.%s.sto'%(prefix,db)
    tmp = "%s.%s.out"%(prefix,db)
    tbl = "%s.%s.tbl"%(prefix,db)

    cmd = "qhmmsearch --cpu %d -A %s -o %s -E 10 --incE 1e-3 "\
          "--tblout %s %s %s"%(ncpu,out,tmp,tbl,hmm,database[db])

    time = datetime.datetime.now().time()

    ## sbatch
    sbatch = "sbatch -N 1 -n %d -o qhmmer_%s.out --wait --wrap="%(ncpu,db)
    sbatch_cmd = "%s\"%s\""%(sbatch,cmd)
    print('\n',"-%s Performing qhmmsearch on %s: "%(time,db),cmd,'\n')
    if not skip_search:
        #subprocess.check_output(cmd,shell=True)
        subprocess.check_output(sbatch_cmd,shell=True)
    else:
        time = datetime.datetime.now().time()
        print('\n',"-%s Skipping qhmmsearch on %s: "%(time,db),'\n')

    ## specify output file names ##
    a3mout = '%s.%s.a3m'%(prefix,db)

    ## if no alignments returned ##
    if is_empty_file(out):
        shutil.copy(a3m,a3mout)
        return 0
    else:
        return 1

def build_and_search_custom_db(prefix,previous_db,db):
    ## fetch sequence from hmmer search table ##
    tbl   = '%s.%s.tbl'%(prefix,db)
    fseqs = '%s.%s.fseqs'%(prefix,db)

    time = datetime.datetime.now().time()
    cmd = "esl-sfetch -f %s %s | sed 's/*//g' > %s"%(database[db],tbl,fseqs)
    sys.stdout.write("-%s Fetch match sequences "\
                     "from hmmer search on %s: %s\n"%(time,db,cmd))

    #if not db=='uniref':
    #if not db=="metaclust":
    #    subprocess.check_output(cmd,shell=True)
    subprocess.check_output(cmd,shell=True)

    ## in case that <fseqs> is empty, copy some content into it
    randseq1 = {"randseq1":"MDINIGTQGGQEHVKMFMYPRNPRVQFRTQVSLMQYSTWTADWYPGA"}
    randseq2 = {"randseq2":"MDAAAGTQGGQKAVDPPPYPRNPRVQFSPPPPDLQYSTAAADWYPGA"}
    randomseq = {**randseq1,**randseq2}
    seq2fasta(randomseq,fseqs)
  
    ## then cluster and build custom hhblits database ##
    fseqs_db = "%s.%s_fseqs_db"%(prefix,db) ## create a sequence database ##
    clust_db = "%s.%s_clust_db"%(prefix,db)
    clust_msa = "%s.%s_clustMsa"%(prefix,db)
    cmd = "mmseqs createdb %s %s"%(fseqs,fseqs_db)

    time = datetime.datetime.now().time()
    sys.stdout.write("-%s Cluster match sequences "
                     "from hmmer search on %s: %s\n"%(time,db,cmd))
    if os.path.isdir("./tmp/"):
        shutil.rmtree("./tmp/")
    subprocess.check_output(cmd,shell=True) #create sequence database#
    
    cmd = "mmseqs createindex %s tmp"%fseqs_db ## create index command
    time = datetime.datetime.now().time()
    sys.stdout.write("-%s Cluster match sequences "\
                     "from hmmer search on %s: %s\n"%(time,db,cmd))
    if os.path.isdir("./tmp/"):
        shutil.rmtree("./tmp/")
    subprocess.check_output(cmd,shell=True) ## create index for the database
    
    cmd = "mmseqs cluster %s %s tmp --min-seq-id 0.3 "\
          "--threads %d"%(fseqs_db,clust_db,ncpu)
    time = datetime.datetime.now().time()
    sys.stdout.write("-%s Cluster match sequences "\
                     "from hmmer search on %s: %s\n"%(time,db,cmd))
    #if not db=="uniref":
    #    subprocess.check_output(cmd,shell=True) ## cluster the sequences
    if os.path.isdir("./tmp/"):
        shutil.rmtree("./tmp/")
    subprocess.check_output(cmd,shell=True) ## cluster the sequences
    
    fseqs_db = "%s.%s_fseqs_db"%(prefix,db) ## create a sequence database ##
    clust_db = "%s.%s_clust_db"%(prefix,db)
    clust_msa = "%s.%s_clustMsa"%(prefix,db)
    cmd = "mmseqs result2msa %s %s %s %s --compressed 1 "\
          "--diff 9999 --max-seq-id 0.99"%(fseqs_db,fseqs_db,clust_db,clust_msa)
    time = datetime.datetime.now().time()
    print(time,cmd)
    if os.path.isdir("./tmp/"):
        shutil.rmtree("./tmp/")
    subprocess.check_output(cmd,shell=True) ## result to MSA
  
  
    cmd ="sbatch -N 1 -n %d -o hhblits_db_translate_%s.out "\
      "--wait --wrap=\"mpirun -np %d cstranslate_mpi "\
      "-i %s -o %s_cs219 "\
      "-f -x 0.3 -c 4 -I ca3m -b\""%(ncpu,db,ncpu,clust_msa,clust_msa)
    time = datetime.datetime.now().time()
    print(cmd)
    sys.stdout.write("-%s Create custom hhblits database  "\
                     "from hmmer search on %s: %s\n"%(time,db,cmd))
    try:
        subprocess.check_output(cmd,shell=True) ## build hhblits database
      
        ## hblits on the clustom database ##
        database[clust_msa] = clust_msa
        ina3m  = "%s.%s.a3m"%(prefix,previous_db)
        a3mout = "%s.a3m"%(clust_msa)
        hhmout = "%s.a3m"%(clust_msa)
      
        time = datetime.datetime.now().time()
        sys.stdout.write("-%s HHBLits search against customDB %s: %s\n"%(time,clust_msa,cmd))
        hhblits(ina3m,a3mout,hhmout,ncpu,clust_msa)
    except subprocess.CalledProcessError:
        ## if hmmer search returns zero results
        ina3m  = "%s.%s.a3m"%(prefix,previous_db)
        a3mout = "%s.a3m"%(clust_msa)
        shutil.copy(ina3m,a3mout)
        

    ## if uniref stage, replace previous if Nf is greater, else use previous ##
    #clust_msa = "%s.%s_clustMsa"%(prefix,db)
    #a3mout = "%s.a3m"%(clust_msa)
    if db=="uniref":
        #chop_seq(a3mout)
        a3m_prev = "%s.%s.a3m"%(prefix,previous_db)

        nf_curr = check_Nf(a3mout)
        nf_prev = check_Nf(a3m_prev)

        a3mout_now = '%s.%s.a3m'%(prefix,db)

        print('\n',"- Previous Nf is %f, now Nf is %f "%(nf_prev,nf_curr),'\n')
        if nf_curr > nf_prev:
            shutil.copy(a3mout,a3mout_now) 
        else:
            shutil.copy(a3m_prev,a3mout_now)

    ## if other database, combine previous a3m and the resulting a3m ##
    elif db=="metaclust" or db=="tara":
        #chop_seq(a3mout)
        a3m_prev = "%s.%s.a3m"%(prefix,previous_db)

        a3m_comb = "%s.%s.combined.a3m"%(prefix,db)
        cmd = "hhalign -id 99 -diff inf -cov 50 "\
              "-i %s -t %s -oa3m %s"%(a3m_prev,a3mout,a3m_comb)
        sys.stdout.write('\n- Combining previous search results: %s\n'%cmd)
        subprocess.check_output(cmd,shell=True)

        ## filter out redundant sequences ##
        a3mout = '%s.%s.a3m'%(prefix,db)
        cmd = "hhfilter -id 99 -diff inf -cov 50 -i %s -o %s"%(a3m_comb,a3mout)
        sys.stdout.write('\n- Filtering redundant sequences: %s\n'%cmd)
        subprocess.check_output(cmd,shell=True)

        ## realign ##
        a3m_realn = '%s.%s.realn.a3m'%(prefix,db)
        a3m_bak   = '%s.%s.a3m.bak'%(prefix,db)
        mafft_realign(a3mout,a3m_realn,ncpu)
        shutil.copy(a3mout,a3m_bak)
        shutil.copy(a3m_realn,a3mout)


def jackhmmer_on_db(prefix,previous_db,db,skip_search=False):
    if skip_search:
        sys.stdout.write("- Skipping search on %s\n"%db)
    else:
        #if not db=='uniref':
        #    status = jackhmmer_search(prefix,previous_db,db,skip_search)
        #else:
        #    status = 1
        status = jackhmmer_search(prefix,previous_db,db,skip_search)
        if status:
            build_and_search_custom_db(prefix,previous_db,db)
        else:
            sys.stdout.write("\n- No results hmmer on %s, Skipping...\n"%db)

def hmmer_on_db(prefix,previous_db,db,skip_search=False):
    if skip_search:
        sys.stdout.write("- Skipping search on %s\n"%db)
    else:
        #if not db=='metaclust':
        #    status = hmmer_search(prefix,previous_db,db,skip_search)
        #else:
        #    status = 1
        status = hmmer_search(prefix,previous_db,db,skip_search)
        if status:
            build_and_search_custom_db(prefix,previous_db,db)
        else:
            sys.stdout.write("\n- No results hmmer on %s, Skipping...\n"%db)

class MSABuilder():
    def __init__(self,query,ncpu=4):
        self.query = query
        self.ncpu  = ncpu
        ## get prefix and other staffs from input filename ##
        prefix     = os.path.splitext(os.path.basename(query))[0]
        self.a3m   = "%s.a3m"%prefix
        self.prefix= prefix

        ## copy query to current dir if its in some other folder ##
        print(self.query)
        query_seq = fasta2seq(self.query)
        key       = list(query_seq)[0]
        seq       = query_seq[key]
        seq       = seq.replace('X','A')
        seq       = seq.replace('U','C')
        query_seq[key] = seq
        query = "%s.fasta"%prefix
        seq2fasta(query_seq,query)

    def search_uniclust30(self,db,skip_search):
        hhblits_on_db(self.prefix,self.ncpu,db,skip_search=skip_search)

    def enough_msa(self,db):
        time = datetime.datetime.now().time()
        sys.stdout.write("-%s Checking if result msa "\
                         "from %s has enough Nf\n"%(time,db))
        return enough_msa(self.prefix,db)

    def build(self,db1,db2,db3,db4,remove_intermediate=True):
        self.search_uniclust30(db1,skip_search=False)
        final_db = db1
        if not self.enough_msa(db1):
            jackhmmer_on_db(self.prefix,db1,db2,skip_search=False)
            final_db = db2
            if not self.enough_msa(db2):
                hmmer_on_db(self.prefix,db2,db3,skip_search=False)
                final_db = db3
                if not self.enough_msa(db3):
                    hmmer_on_db(self.prefix,db3,db4,skip_search=False)
                    final_db = db4

        a3mout = '%s.%s.a3m'%(self.prefix,final_db)
        shutil.copy(a3mout,self.a3m)


        if os.path.isdir("./tmp/"):
            shutil.rmtree("./tmp/")
        if remove_intermediate==True:
            for f in glob.glob("%s*"%self.prefix):
                dbs = [db1,db2,db3,db4]
                if any([s in f for s in dbs]) or f.endswith(".hhr"):
                    os.remove(f)


if __name__ == "__main__":

    #query = "seq020.fa"
    query = sys.argv[1]
    ncpu  = 4 

    db1 = "uniclust"
    db2 = "uniref"
    db3 = "tara"
    db4 = "metaclust"
    
    msa = MSABuilder(query,ncpu=ncpu)
    msa.build(db1,db2,db3,db4,remove_intermediate=True) 
