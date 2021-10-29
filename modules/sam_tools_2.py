import pandas as pd
import os
import re
import pysam
import copy
import glob

class junc_site():
    def __init__(self) -> None:
        self.read5_chr=""
        self.read5_start=0
        self.read5_end=0
        self.read3_chr=""
        self.read3_start= 0
        self.read3_end = 0
        # the other end of fragaments
        self.reverse_dic={
            "f":0,
            "r":0
        }
        
        self.fragment_num = 0
        # it is better use birary treee to sort
        self.readx_ls = []
        self.readx_loc_dic={
            "in":0,
            "out":0
        }

    # (read5_start)>>>>>>>(read5_end)|(read3_start)>>>>>>>>(read3_end) ------readx----------
    def getSiteID(self):
        return "{}_{}!{}_{}".format(self.read5_chr,self.read5_end,self.read3_chr, self.read3_start)

    def getLength(self):
        if self.read3_chr==self.read5_chr:
            return abs(self.read3_start-self.read5_end)+1
        else:
            return -1
    def getReadXlocs(self):
        if self.read5_chr!=self.read3_chr:
            return [False]*len(self.readx_ls)
        result=[]
        for readx in self.readx_ls:
            readx_start = readx["start"]
            readx_end = readx["end"]
            
            if readx_start>readx_end: readx_start,readx_end = readx_end,readx_start
            if self.read5_end>self.read3_start:
                if readx_start>self.read3_start and readx_end<self.read5_end:
                    result.append(True)
                else:
                    result.append(False)
            else: #self.read3_start>self.read5_end
                if readx_start>self.read5_end and readx_end<self.read3_start:
                    result.append(True)
                else:
                    result.append(False)
        return result

    
    # self validate
    def validate(self):
        return self.read3_chr!="" and self.read5_chr!=""

    @staticmethod
    def toCSVTitle():
        return ",".join(["id","read5_chr","read5_end","read3_chr","read3_start","len","frag_num","f_frag_num","r_frag_num","readx_in","readx_out"])
        
    def toCSV(self):
        return ("{},{},{},{},{},{},{},{},{},{},{}".format(self.getSiteID(),\
                                                self.read5_chr,self.read5_end,\
                                                self.read3_chr,self.read3_start,\
                                                self.getLength(), self.fragment_num,\
                                                self.reverse_dic["f"],self.reverse_dic["r"],\
                                                self.readx_loc_dic["in"],self.readx_loc_dic["out"]\
                                                ))
    def toBED(self):
        if (self.read3_chr!=self.read5_chr):
            return ""
        else:
            return "{}\t{}\t{}\t{}\t{}\t+".format()

#bam_file="/home/hqyone/mnt/2tb/eccDNA/data-final/eccdna/test_sample_id/bam/test_sample_id_hg38_sort_du_byname.bam"
bam_file="/home/hqyone/mnt/2tb/eccDNA/data-final/eccdna/D21010704/bam/D21010704_hg38_sort_du_byname.bam"
js_csv = "/home/hqyone/mnt/2tb/eccDNA/data-final/eccdna/test_sample_id/bam/test_sample_id_hg38_sort_du_byname.csv"


def getJunCSV(bam_file, js_csv):
    if not os.path.isfile(bam_file+".bai"):
        fp  = open(bam_file+".bai",'x')
        fp.close()
    BAM = pysam.AlignmentFile(bam_file, 'rb')
    q_name =""
    read_ls = []
    js_dic = {}
    sam_total_frag_num=0
    for read in BAM:
        if (q_name!="" and q_name!=read.query_name):
            sam_total_frag_num+=1
            if (len(read_ls)>2):
                js = junc_site()
                readx_is_reverse = None
                find_readx=False
                for r in read_ls:
                    if (re.match(r"^(\d+)M$",r.cigarstring)):
                        js.readx_ls.append({
                            "chr":r.reference_name,
                            "start":r.reference_start,
                            "end":r.reference_end
                        })
                        readx_is_reverse = r.is_reverse
                        find_readx = True
                        break
                if find_readx:
                    for r in read_ls:
                        chrom=r.reference_name
                        start=r.reference_start # In bam file the start pos is always left residue
                        end=r.reference_end
                        seq = r.seq
                        rev=r.is_reverse
                        cigar_str = r.cigarstring
                        #print(read.query_name+"---------"+r.cigarstring)
                        #print("{}:{}-{}--{}--{}".format(chrom,start, end, seq, rev))
                        if rev!=readx_is_reverse:
                            m=re.match(r"^(\d+)[S](\d+)M$",cigar_str) #5'(+, -)
                            n=re.match(r"^(\d+)M(\d+)[S]$",cigar_str) #3'(-, +)
                            if m:
                                match_len=int(m.groups(0)[1])
                                #print("{}:{}:{}".format(start,match_len,rev))
                                js.read3_start = start
                                js.read3_chr = chrom
                                js.read3_end = end 
                            elif n:
                                match_len=int(n.groups(0)[0])
                                #print("{}:{}:{}".format(end,match_len,rev))
                                js.read5_start = start
                                js.read5_chr = chrom
                                js.read5_end = end
                    if (js.validate()):
                        js_id = js.getSiteID()
                        readx_loc_in=js.getReadXlocs()[0]
                        if js_id not in js_dic:
                            js.fragment_num=1
                            js_dic[js_id] = js
                        else:
                            js_dic[js_id].fragment_num+=1
                        if readx_is_reverse:
                            js_dic[js_id].reverse_dic["f"]+=1
                        else:
                            js_dic[js_id].reverse_dic["r"]+=1
                        if readx_loc_in:
                            js_dic[js_id].readx_loc_dic["in"]+=1
                        else:
                            js_dic[js_id].readx_loc_dic["out"]+=1
                    #print("{}:{}:{}".format(js_id, js.getLength(),js_dic[js_id].fragment_num))

                if len(js_dic.keys())>2000:
                    break
            read_ls=[]
        read_ls.append(read)
        q_name = read.query_name
    print("fragment total num: {}\n".format(sam_total_frag_num))
    with open(js_csv,'w') as CSV:
        CSV.write(js.toCSVTitle()+"\n")
        for key, js in js_dic.items():
            CSV.write(js.toCSV()+'\n')

# find all bam files in the data-final folder
sample_sheet = '/home/hqyone/mnt/2tb/eccDNA/data-raw/eccdna_sample_sheet.csv'
s_df = pd.read_csv(sample_sheet)
out_dir = '/home/hqyone/mnt/2tb/eccDNA/data-final/eccdna/jscsv'
data_final_dir = '/home/hqyone/mnt/2tb/eccDNA/data-final/'
bam_ls = glob.glob('{}/eccdna/**/*_byname.bam'.format(data_final_dir), recursive=True)
for  bam in bam_ls:
    if "test" in bam:
        continue
    s_id = bam.split("/")[-3]
    s_name = list(s_df.loc[s_df["sample_id"]==s_id]["sample_name"])[0].replace('_eccdna',"")
    js_csv = "{}/{}_{}_js.csv".format(out_dir,s_id,s_name)
    print(s_id,s_name, js_csv)
    getJunCSV(bam_file, js_csv)
