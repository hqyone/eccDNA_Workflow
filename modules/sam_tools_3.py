import pandas as pd
import os
import re
import pysam
import copy
import glob


class junc_site():
    def __init__(self) -> None:
        self.read5=None
        self.read3=None
        self.readx=None
        # Statisic field
        self.frag_num=0
        self.closed_frag_num=0
        self.read5_rev_num=0
        self.read3_rev_num=0
        
            
        

    # (read5_start)>>>>>>>(read5_end)|(read3_start)>>>>>>>>(read3_end) ------readx----------
    def getSiteID(self):
        return f"{self.read5.reference_name}_{self.read5.reference_end}!"+\
               f"{self.read3.reference_name}_{self.read3.reference_start}"
        #read5_isrev = 1 if self.read5.is_reverse else 0
        read3_isrev = 1 if self.read3.is_reverse else 0

    @staticmethod
    def compareRead(read1, read2):
        '''
            return True if read1 located later than read2 else False
        '''
        if read1.reference_name==read2.reference_name:
            return read1.reference_start > read2.reference_start
        else:
            return read1.reference_name > read2.reference_name

    def isClosed(self):
        if self.read5 and self.read3:
            if (self.read5.reference_name  == self.read3.reference_name) and (self.read5.is_reverse == self.read3.is_reverse):
                if self.read5.reference_end > self.read3.reference_start:
                        return True
        return False
    
    # self validate and initialization
    def initilized(self):
        if self.read5 and self.read3:
            if self.read3.reference_name!="" and self.read5.reference_name!="":
                self.frag_num=1
                if self.isClosed(): self.closed_frag_num=1
                if self.read5.is_reverse: self.read5_rev_num=1
                if self.read3.is_reverse: self.read3_rev_num=1
                return True
        return False

    def merge(self, js):
        if self.getSiteID() == js.getSiteID():
            self.frag_num+=1
            self.closed_frag_num +=js.closed_frag_num
            self.read5_rev_num +=js.read5_rev_num
            self.read3_rev_num +=js.read3_rev_num
    
    def getLength(self):
        if self.read5 and self.read3:
            if self.read5.reference_name==self.read3.reference_name:
                return abs(self.read3.reference_start-self.read5.reference_end)+1
        return -1
    

    @staticmethod
    def toCSVTitle():
        return ",".join(["id","read5_chr","read5_end","read3_chr","read3_start","len","frag_num","close_frg_num","read5_rev_num","read3_rev_num"])
        
    def toCSV(self):
        return (f"{self.getSiteID()},{self.read5.reference_name},{self.read5.reference_end},"+\
                f"{self.read3.reference_name},{self.read3.reference_start},{self.getLength()},{self.frag_num},"+\
                    f"{self.closed_frag_num},{self.read5_rev_num},{self.read3_rev_num}")
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
                readx_is_read1=None
                find_readx=False
                for r in read_ls:
                    if (re.match(r"^(\d+)M$",r.cigarstring)):
                        js.readx=r
                        readx_is_reverse = r.is_reverse
                        readx_is_read1 = r.is_read1
                        find_readx = True
                        break
                if find_readx:
                    for r in read_ls:
                        cigar_str = r.cigarstring
                        if r.is_read1!=readx_is_read1:  # r and readx not from the sample fastq file (tow terminates of a fragment)
                            m=re.match(r"^(\d+)[S](\d+)M$",cigar_str) #5'(+, -)
                            n=re.match(r"^(\d+)M(\d+)[S]$",cigar_str) #3'(-, +)
                            if m:
                                match_len=int(m.groups(0)[1])
                                js.read3 = r
                            elif n:
                                match_len=int(n.groups(0)[0])
                                #print("{}:{}:{}".format(end,match_len,rev))
                                js.read5 = r
                    if (js.initilized()):
                        js_id = js.getSiteID()
                        if js_id not in js_dic:
                            js_dic[js_id] = js
                        else:
                            js_dic[js_id].merge(js)
                    #print("{}:{}:{}".format(js_id, js.getLength(),js_dic[js_id].fragment_num))
                # if len(js_dic.keys())>2000:
                #     break
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
    getJunCSV(bam, js_csv)
