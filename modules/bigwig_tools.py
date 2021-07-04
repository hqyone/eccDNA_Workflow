'''
Author: Quanyuan(Leo) He
Email: hqyone@gmail.com
Insititute: Hunan Normal Univeristy
Date: 1969-12-31 17:00:00
LastEditTime: 2021-04-17 22:08:06
LastEditors: Quanyuan(Leo) He
Description: 
FilePath: /eccDNA/code/source/bigwig_tools.py
License: The MIT License (MIT)
'''
import pyBigWig
#from bitarray import bitarray



def refineBedFile(bed, bigwig, outbed, offset=400, cutoff=30):
    bw = pyBigWig.open(bigwig)
    with open(bed, 'r') as INPUT, open(outbed, 'w') as OUTPUT:
        for line in  INPUT:
            contents = line.strip().split("\t")
            chr= contents[0]
            start = int(contents[1])
            end = int(contents[2])
            remain_str = "\t".join(contents[3:])
            (out_start, out_end) = updateRegion(chr, start, end, bw, offset, cutoff)
            OUTPUT.write(chr+"\t"+str(out_start)+"\t"+str(out_end)+"\t"+remain_str+"\n")
                
def updateRegion(chr, start, end, bw, offset=400, cutoff=30):
    out_start = start
    out_end = end
    if (chr in bw.chroms() and chr!="chrM"):
        left_values = bw.values(chr, start-offset, start)
        for i in range(offset-1):
            if (left_values[-(i+1)]>cutoff):
                out_start-=1
            else:
                break
        
        right_values = bw.values(chr, end, end+offset)
        for i in range(offset-1):
            if (right_values[i]>cutoff):
                out_end+=1
            else:
                break
    return (out_start, out_end)

def findRegin(chr, start, direct, bw, offset=10000, cutoff=30):
    out_start = start
    out_end = start
    if (chr in bw.chroms()):
        if direct=="+":
            values = bw.values(chr, start, start+offset)
            for i in range(offset-1):
                if (values[i]>cutoff):
                    out_end+=1
                else:
                    break
        else:
            values = bw.values(chr, start-offset, start)
            for i in range(offset-1):
                if (values[-(i+1)]>cutoff):
                    out_start-=1
                else:
                    break
    return (out_start, out_end)

class JoinSiteRegion:
    def __init__(self, chr, start, end, direction, for_reads, rev_reads, coverage):
        self.chr=chr
        self.start=start
        self.end=end
        self.direction = direction
        self.coverage=coverage
        self.for_reads= for_reads
        self.rev_reads = rev_reads
    
    def getID(self):
        return self.chr+":"+str(self.start)+"-"+str(self.end)+":"+str(self.direction)
    
    def isValidated(self):
        if (abs(self.start-self.end)<=2):
            return False
        else:
            return True
        

def findRegion2(bw_array, chr, start, direct, for_reads, rev_reads, cutoff=30):
    out_start = start
    out_end = start
    #region_str=""
    coverage = 0
    sum = 0
    if direct=="+":
        while out_end<len(bw_array):
            if bw_array[out_end-1]>cutoff:
                sum +=bw_array[out_end-1]
                out_end+=1
            else:
                break;
    else:
        while out_start-1>=0:
            if bw_array[out_start-1]>cutoff:
                sum+=bw_array[out_start-1]
                out_start-=1
            else:
                break
    if (out_start!=out_end):
        coverage = int(sum/(abs(out_end - out_start)))
    #if (out_start!=out_end):
    #    region_str = str(out_start)+"-"+str(out_end)+"-"+direct
    jsr = JoinSiteRegion(chr, out_start, out_end, direct, for_reads, rev_reads, coverage)
    return jsr

def bw2dic(bwfile, genomeSize_file, cutoff=30):
    seq_dic={}
    bw = pyBigWig.open(bwfile)
    with open(genomeSize_file, 'r') as SIZE:
        for l in SIZE:
            contents = l.strip().split("\t")
            chr = contents[0]
            size = int(contents[1])
            if "_" not in chr:
                temp_ls=bw.values(chr, 0, size-1)
                #bit_ls = bitarray(len(temp_ls))
                bit_ls =[False]*len(temp_ls)
                for i in range(0,len(temp_ls)):
                    if temp_ls[i]>=cutoff:
                      bit_ls[i] = True
                temp_ls=None
                seq_dic[chr]=bit_ls
    return seq_dic

# Test
# bed_file="/home/hqyone/mnt/2t/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak"
bigwig_file= "/media/hqyone/2tb1/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs/sample1.bigwig"
# refined_bed_file="/home/hqyone/mnt/2t/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.refine"

# bw = pyBigWig.open(bigwig_file)
# a =bw.values("chr12",104327580,104327600)
# b =bw.values("chr12",104327205,104327605)
# updateRegion("chr12",104327605, 104329387, bw)

#refineBedFile(bed_file, bigwig_file, refined_bed_file)

genome_size_file = "/media/hqyone/2tb1/eccDNA/genome/chrom_size/hg38.chrom.sizes"
#dic =  bw2dic(bigwig_file, genome_size_file)
print("hqyone")
