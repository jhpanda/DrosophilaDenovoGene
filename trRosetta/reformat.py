
import sys
from fasta2seq import *

def extract_regions(pref,start,end):
    a3m = '%s.a3m'%pref
    fas = '%s_trunc.a3m'%pref
    sequence = fasta2seq(a3m)

    flag = 0
    newsequence = {}
    for key in sequence:
        seq = sequence[key]
        if start!=None and end!=None:
            newseq = seq[start-1:end]
        else:
            newseq = seq
        print(key,seq)
        newsequence[key] = newseq
        flag += 1

    f = open(fas,"w")
    for key in newsequence:
        f.write("%s\n"%newsequence[key])
    f.close()
    #seq2fasta(newsequence,fas,linebreak=False)

def main(pref,start,end):
    extract_regions(pref,start,end)
    #if start!=None and end!=None:
    #    extract_regions(pref,start,end)
    
#
if __name__ == '__main__':
    if len(sys.argv)<2:
        print("usage: python reformat.py pref start(optional) end(optional)")
        sys.exit(0)

    pref = sys.argv[1]
    if len(sys.argv)>2:
        start = int(sys.argv[2])
        end   = int(sys.argv[3])
    else:
        start = None
        end   = None
    main(pref,start,end)
