#!/usr/bin/python3

# search for boundary using bigwig and macs bed file
import os, sys, re
from typing import Sequence
import pysam
from pysam import AlignedSegment


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
def rc(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))

class JuncSite2():
    def __init__(self):
        self.l_chr =""
        self.l_loc=0
        self.l_direct = "+"
        self.l_seq=""
        self.r_chr = ""
        self.r_loc=0
        self.r_direct = "+"
        self.r_seq= ""
        self.for_read_num=0
        self.rev_read_num=0
        
    def getID(self):
        return self.l_chr+"_"+str(self.l_loc)+"_"+self.l_direct+"#"+self.r_chr+"_"+str(self.r_loc)+"_"+self.r_direct

    def CombineSite(self, j2):
        if self.getID()==j2.getID():
            self.for_read_num+=j2.for_read_num
            self.rev_read_num+=j2.rev_read_num        
        
        
# The juncSite is abstract object to descript junction region
class JuncSite():
    def __init__(self, ji):
        self.id = ji.getID()
        self.type=ji.type
        self.chr=ji.chr
        self.loc=ji.loc #always left hand
        self.seq_num=1
        self.left_seq_dic={}
        self.left_seq = ""
        self.left_coverage = 0
        self.right_seq_dic={}
        self.right_seq = ""
        self.right_coverage = 0
        self.left_next_map_locs={}
        self.right_next_map_locs={}
        self.rev_read_num = 0
        self.for_read_num = 0
        self.appendRead(ji)
    
    def getSupportReadNumber(self):
        left_support_reads_num = 0
        right_support_reads_num=0
        for seq, item in self.left_seq_dic.items():
            left_support_reads_num+=item
        for seq, item in self.right_seq_dic.items():
            right_support_reads_num+=item
        total = left_support_reads_num+right_support_reads_num
        return (total, left_support_reads_num, right_support_reads_num)   
            
    
    def toBedStr(self):
        (t, l, r) = self.getSupportReadNumber()
        return f"{self.chr}\t{self.loc}\t{self.loc+1}\t{self.id}\t{t}\t{l}\t{r}"
        

    def appendRead(self,ji):
        if (ji.getID()==self.id):
            self.seq_num+=1
            if ji.left_seq not in self.left_seq_dic:
                self.left_seq_dic[ji.left_seq]=0
            self.left_seq_dic[ji.left_seq]+=1
            if ji.right_seq not in self.right_seq_dic:
                self.right_seq_dic[ji.right_seq]=0
            self.right_seq_dic[ji.right_seq]+=1
            if ji.left_next_map_loc!="":
                if ji.left_next_map_loc not in self.left_next_map_locs:
                    self.left_next_map_locs[ji.left_next_map_loc]=0
                self.left_next_map_locs[ji.left_next_map_loc]+=1
            if ji.right_next_map_loc!="":
                if ji.right_next_map_loc not in self.right_next_map_locs:
                    self.right_next_map_locs[ji.right_next_map_loc]=0
                self.right_next_map_locs[ji.right_next_map_loc]+=1
            if ji.rev:
                self.rev_read_num+=1
            else:
                self.for_read_num+=1
        

    def getSiteConsenseSeq(self, min_ratio=0.4, min_reads=3):
        l_Seq_ls=self.left_seq_dic.keys()
        r_Seq_ls=self.right_seq_dic.keys()
        (self.left_seq, self.left_coverage) = self.getConsenseSeq(l_Seq_ls, "right", min_ratio, min_reads)
        (self.right_seq, self.right_coverage) = self.getConsenseSeq(r_Seq_ls,"left",min_ratio, min_reads)

    def getConsenseSeq(self, seq_ls, aligned_to, min_ratio=0.5, min_reads=15):
        cons_seq=""
        cons_deep=0
        for i in range(200):
            cur_char_dic={}
            read_sum=0
            max_num=0
            max_char=""
            for seq in seq_ls:
                if i<len(seq):
                    char = seq[-1*(i+1)] 
                    if (aligned_to == "left"): char=seq[i]
                    if char not in cur_char_dic:
                        cur_char_dic[char]=0
                    read_sum+=1
                    cur_char_dic[char]+=1
                    if cur_char_dic[char]>max_num:
                        max_num=cur_char_dic[char]
                        max_char=char
            if (read_sum==0): break
            max_ratio= max_num/float(read_sum)
            if max_ratio<min_ratio or max_num<min_reads:break
            else:
                if (aligned_to == "left"):
                    cons_seq=cons_seq+max_char
                else:
                    cons_seq=max_char+cons_seq
                if cons_deep<max_num: cons_deep=max_num
        return (cons_seq, cons_deep)

    def testConsenseSeq(self):
        seq_ls=["ATTTAGGGGGTASSTSSKT","AATTAGGGGGTASSTSS","AATTAGGGGGTAS","AATTAGGGGG"]
        print(self.getConsenseSeq(seq_ls, "left", 0.4, 2))
        seq_ls=["ATTTAGGGGGTASSTSSKT","AATTAGGGGGTASSTSSKT","AATTAGGGGGGTASSTSSKT","AATTAGGGGG"]
        print(self.getConsenseSeq(seq_ls, "right", 0.4, 2))                  

class JuncInfor():
    def __init__(self):
        self.type=""
        self.chr=""
        self.loc="" #always left hand and means the cliped sites
        self.next_map_loc=""
        self.left_seq="" # All transform to plus strand
        self.right_seq="" # All transform to plus strand
        self.left_next_map_loc="" #
        self.right_next_map_loc="" #
        self.rev = False
    
    def toStr(self):
        return f"{self.chr},{self.type},{self.loc},{self.left_seq},{self.right_seq},{self.next_map_loc}"
    
    def getID(self):
        return f"{self.chr}_{self.loc}_{self.type}"

def getJunctionInfor(read):
    if (read==None):
        return None
    cigar_str=read.cigarstring
    chrom=read.reference_name
    start=read.reference_start # In bam file the start is left residue
    end=read.reference_end
    query_seq=read.query_sequence
    rev=read.is_reverse
    if not cigar_str: return None
    if (chrom=="chr1" and end==14953):
        print("hqyone")
    m=re.match(r"^(\d+)[S](\d+)M$",cigar_str) #5'(+, -)
    n=re.match(r"^(\d+)M(\d+)[S]$",cigar_str) #3'(-, +)
    if m or n:
        ji=JuncInfor()
        ji.chr=chrom
        if m: 
            left_len=int(m.groups(0)[0])
            ji.loc=start
            ji.type=5
            ji.seq=query_seq
            ji.left_seq=query_seq[0:left_len]
            ji.right_seq=query_seq[left_len:]
            ji.rev = rev
        elif n:
            left_len=int(n.groups(0)[0])
            ji.loc=start+left_len
            ji.type=3
            ji.seq = query_seq
            ji.left_seq=query_seq[0:left_len]
            ji.right_seq=query_seq[left_len:]
            ji.rev = rev
        # strand="+"
        # if read.is_reverse: strand='-'
        # if m and strand=="+":
        #     left_len=int(m.groups(0)[0])
        #     ji.loc=start-1
        #     ji.type=5
        #     ji.left_seq=query_seq[0:left_len]
        #     ji.right_seq=query_seq[left_len:]
        # elif m and strand=="-":
        #     left_len=int(m.groups(0)[0])
        #     ji.loc=end
        #     ji.type=3
        #     ji.left_seq=rc(query_seq[left_len:])
        #     ji.right_seq=rc(query_seq[0:left_len])
        # elif n and strand=="+":
        #     left_len=int(n.groups(0)[0])
        #     ji.loc=end
        #     ji.type=3
        #     ji.left_seq=query_seq[0:left_len]
        #     ji.right_seq=query_seq[left_len:]
        # elif n and strand=="-":
        #     left_len=int(n.groups(0)[0])
        #     ji.loc=start-1
        #     ji.type=5
        #     ji.left_seq=rc(query_seq[left_len:])
        #     ji.right_seq=rc(query_seq[0:left_len:])
        # else:
	    #     return None
        return ji
    else:
        return None

# Mimic pysam read for testing
class myAlignmentSeg():
    def __init__(self, ciga, ref_name,strand,start, end, seq) -> None:
        self.cigarstring=ciga
        self.reference_name = ref_name
        self.reference_start = start
        self.reference_end = end
        self.query_sequence = seq.upper()
        self.strand=strand
        self.is_reverse=(strand=="-")

def test():
    #       3|4    10|11
    #  5'-ATG|AAT.AAG|GAC-3' 
    #  3'-TAC|TTA.TTC|CTG-5'
    read1=myAlignmentSeg("3S3M","chr1","+",4,6,"ATGAAT")
    if (getJunctionInfor(read1)):print(getJunctionInfor(read1).toStr())

    read2=myAlignmentSeg("3M3S","chr1","+",8,10,"AAGGAC")
    if (getJunctionInfor(read2)):print(getJunctionInfor(read2).toStr())

    # chr1,5,3,ATG,AAT
    # chr1,3,10,AAG,GAC

test()
js = JuncSite(JuncInfor())
js.testConsenseSeq()