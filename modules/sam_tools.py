
from modules.junction_infor import JuncInfor, JuncSite, getJunctionInfor,JuncSite2
from modules.blast_tools import RunBLASTN
from modules.bigwig_tools import refineBedFile, findRegin, findRegion2, JoinSiteRegion
from modules.kosaraju_graph import Graph
import pandas as pd
import pyBigWig
import os
import re
import pysam
import copy

test_dir = "/home/hqyone/mnt/2tb"
sample_bam=test_dir+"/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704/bam/sample1_ns_sort_du.bam"
ref_fa = test_dir+"/eccDNA/data-raw/genome/hg38_p13/GRCh38.p13.genome.fa"
out_dir = test_dir+"/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704"

# A fragment may contains two reads （Read1 and Read2)
class fragment():
    def __init__(self, name) -> None:
        self.query_name=name
        self.read1=None
        self.read2=None
        self.lenth=0

    def addRead(self, read):
        if (read.query_name==self.query_name):
            if (read.is_read1):
                self.read1=read
            elif(read.is_read2):
                self.read2=read
            self.length = read.template_length
    
    def isFull(self):
        return self.read1 and self.read2
    
    def createJuncInfor_ls(self):
        ji_ls=[]
        read1_juninfor = getJunctionInfor(self.read1)
        read2_juninfor = getJunctionInfor(self.read2)
        if read1_juninfor and self.read2:
            read1_juninfor.right_next_map_loc=f"{self.read2.reference_name}:{self.read2.reference_start}-{self.read2.reference_end}"
            ji_ls.append(read1_juninfor)
        if read2_juninfor and self.read1:
            read2_juninfor.left_next_map_loc=f"{self.read1.reference_name}:{self.read1.reference_start}-{self.read1.reference_end}"
            ji_ls.append(read2_juninfor)        
        return ji_ls
        
            
    def isSpanRead(self):
        if self.isFull():
            return self.read1.reference_id!=self.read2.reference_id
        else:
            return False

# Read sam file which is sorted by query_name
def getJunSiteDic(sam, test_read_num=-1):
    BAM = pysam.AlignmentFile(sam, 'rb')
    #for read in samfile.fetch('chr1', 100, 120):
    junc_site_dic={}
    index=0
    frag=None
    query_name=""
    for read in BAM:
        index+=1
        if index % 100000==0:
            print(index)
        if (test_read_num>0 and index>test_read_num):
           break
        if query_name=="": 
            frag = fragment(read.query_name)
            frag.addRead(read)
        elif (read.query_name!=query_name):
            # it is new read_pair
            junc_infor_ls = frag.createJuncInfor_ls()
            for ji in junc_infor_ls:
                ji_id = ji.getID()
                if ji_id not in junc_site_dic:
                    junc_site_dic[ji_id]=JuncSite(ji)
                else:
                    junc_site_dic[ji_id].appendRead(ji)
            frag = fragment(read.query_name)
            frag.addRead(read)
        else:
            frag.addRead(read)
        query_name = read.query_name
    junc_infor_ls = frag.createJuncInfor_ls()
    for ji in junc_infor_ls:
        ji_id = ji.getID()
        if ji_id not in junc_site_dic:
            junc_site_dic[ji_id]=JuncSite(ji)
        else:
            junc_site_dic[ji_id].appendRead(ji)
    for ji_id, site in junc_site_dic.items():
        site.getSiteConsenseSeq()
        
    return junc_site_dic

def JunSiteDicToFasta(junc_site_dic, fasta, min_cov=8, min_len=18):
    with open(fasta,"w") as FASTA:
        for id, site in junc_site_dic.items():
            if site.seq_num>=min_cov:
                (cons_seq, coverage) = (site.left_seq, site.left_coverage)
                if site.type==3:
                    (cons_seq, coverage) = (site.right_seq, site.right_coverage)
                if (coverage>=min_cov and len(cons_seq)>=min_len):
                    FASTA.write(">"+id+"#"+str(coverage)+"\n")
                    FASTA.write(cons_seq+"\n")

# Got the best hits when multiple BLAST result exist            
def compareHit(hit, new_hit, site_id):
    chr=site_id.split("_")[0]
    loc = int(site_id.split("_")[1])
    type = int(site_id.split("_")[2])
    if hit["chr"]== new_hit["chr"]:
        if abs(hit["map_len"]-new_hit["map_len"])<=2:
            if (abs(new_hit["chr_start"]-loc))<abs((hit["chr_start"]-loc)):
                return new_hit
    else:
        if hit["chr"]!=chr and abs(hit["map_len"]-new_hit["map_len"])<=2:
            if type==5:
                if (new_hit["site_read_end"]==len(hit["qseq"])):
                    return new_hit
            elif type==3:
                if (new_hit["site_read_start"]==1):
                    return new_hit
            if (new_hit["chr"]==chr):
                return new_hit
            if new_hit["evalue"]<hit["evalue"]:
                return new_hit
    return hit

def getSiteHitPairStr(site_id, hit):
    hit_seq_str = hit["qseq"]+"_"+str(hit["site_read_start"])+"_"+str(hit["site_read_end"])
    hit_loc_str = hit["chr"]+"_"+str(hit["chr_start"])+"_"+str(hit["chr_end"])+"_"+hit["direction"]
    return site_id+":"+hit_loc_str  #+":"+hit_seq_str

# Here got the extension direction
def creatJuncSite2(js, hit):
    js2 = JuncSite2()
    if js.type == 3:
        js2.l_loc = js.loc
        js2.l_chr = js.chr
        js2.l_direct="-"
        js2.l_seq = js.left_seq
        js2.l_coverage = js.left_coverage
        
        js2.r_chr = hit["chr"]
        js2.r_seq = js.right_seq
        js2.r_coverage = js.right_coverage
        js2.r_loc = hit["chr_start"]
        # hit["chr_start"] is join site location
        if hit["chr_end"]>=hit["chr_start"]:
            js2.r_direct="+"
        else:
            js2.r_direct ="-"
    elif js.type == 5:
        js2.r_loc = js.loc
        js2.r_chr = js.chr
        js2.r_direct="+"
        js2.r_seq = js.right_seq
        js2.l_chr = hit["chr"]
        js2.l_loc = hit["chr_end"]
        # hit["chr_end"] is join site location
        if hit["chr_end"]>=hit["chr_start"]:
            js2.l_direct="-"
        else:
            js2.l_direct ="+"
    js2.for_read_num =js.for_read_num 
    js2.rev_read_num =js.rev_read_num 
    return js2
               
    

def ParseBlastResult(output_tab, junc_site_dic):
    site_hit_dic={}
    Site2_dic = {}
    with open(output_tab,'r') as BLAST_OUT:
        # -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send qseq sseq qlen slen evalue"
        cur_site_id=""
        best_hit={}
        for line in BLAST_OUT:
            contents = line.strip().split("\t")
            site_id = contents[0].split("#")[0] #qseqid
            chr = contents[1] #sseqid
            cur_percent = float(contents[2]) #pident
            cur_length = int(contents[3]) #length
            cur_mismatch = int(contents[4])  #mismatch
            cur_gap = int(contents[5]) #gaps
            site_read_start = int(contents[6]) #qstart
            site_read_end = int(contents[7])  #qend
            chr_start = int(contents[8]) #sstart
            chr_end = int(contents[9])  #send
            qseq = contents[10]  # qseq
            sseq = contents[11]  # sseq
            qlen = int(contents[12])  # qlen
            slen = int(contents[13])  # slen
            evalue = float(contents[14])  # evalue
            direction = "+"
            if chr_start>chr_end:
                direction="-"
            
            if chr_start==499835:
                print("hqyone")
            # blast hit
            hit={
                "evalue":evalue,
                "chr":chr,
                "percent":cur_percent,
                "chr_start":chr_start,
                "chr_end":chr_end,
                "map_len":cur_length, #abs(chr_start-chr_end)+1,
                "direction":direction,
                "site_read_start":site_read_start,
                "site_read_end":site_read_end,
                "qseq":qseq
            }
            if cur_site_id=="":
                best_hit=hit
                cur_site_id=site_id
            elif (cur_site_id!=site_id):
                # new map information about 
                site_hit_dic[cur_site_id]=copy.deepcopy(best_hit)
                if cur_site_id in junc_site_dic:
                    js = junc_site_dic[cur_site_id]
                    js2 = creatJuncSite2(js, best_hit)
                    if js2.getID() not in Site2_dic:
                        Site2_dic[js2.getID()] = js2
                    else:
                        Site2_dic[js2.getID()].CombineSite(js2)
                    
                    print('"'+getSiteHitPairStr(cur_site_id,best_hit)+'",')
                    cur_site_id=site_id
                    best_hit=hit
            else:
                best_hit=compareHit(best_hit, hit, site_id)
        site_hit_dic[cur_site_id]=copy.deepcopy(best_hit)
        if cur_site_id in junc_site_dic:
            js = junc_site_dic[cur_site_id]
            js2 = creatJuncSite2(js, best_hit)
            if js2.getID() not in Site2_dic:
                Site2_dic[js2.getID()] = js2
            else:
                Site2_dic[js2.getID()].CombineSite(js2)
            
        print('"'+getSiteHitPairStr(cur_site_id,best_hit)+'",')
    return (site_hit_dic, Site2_dic)

def JuncSiteToBed(site_hit_dic, bedfile):
    bedtools="/home/hqyone/anaconda3/envs/eccDNA/bin/bedtools"
    with open(bedfile,'w') as BED:
        for id, hit in site_hit_dic.items():
            BED.write(hit.toBedStr()+"\n")
    os.system(f"{bedtools} sort -chrThenScoreD -i {bedfile}>{bedfile+'_sort.bed'}")

def BLASTNConsenseSeq(blastn, mkdb, id, db_fasta, cons_fasta, out_dir, hit_number=20):
    return RunBLASTN(blastn, mkdb, id, db_fasta, cons_fasta, out_dir, eval=0.001, idcutoff=98, hit_number=hit_number)    


def Bed2Dic(bed):
    dic = {}
    with open(bed, "r") as BED:
        for line in BED:
            contents = line.split('\t')
            chr = contents[0]
            start = int(contents[1])
            end = int(contents[2])
            if end<start:
                (end, start)=(start, end)
            if chr not in dic:
                dic[chr]=[]
            dic[chr].append({"s":start, "e":end}) 
    return dic       
        
def annoSiteByBed(site_dic, bed, outfile, offset=50):
    dic = Bed2Dic(bed)
    with open(outfile, 'w') as OUTPUT:
        for i,s in site_dic.items():
            id = s.getID()
            l_region=""
            r_region=""
            if s.l_chr in dic:
                for i in dic[s.l_chr]:
                    if (s.l_loc >i["s"]-offset and s.l_loc<i["e"]+offset):
                        l_region= s.l_chr+":"+str(i["s"])+"-"+str(i["e"])
            if s.r_chr in dic:
                for i in dic[s.r_chr]:
                    if (s.r_loc >i["s"]-offset and s.r_loc<i["e"]+offset):
                        r_region= s.r_chr+":"+str(i["s"])+"-"+str(i["e"])
            OUTPUT.write(id+"\t"+l_region+"\t"+r_region+"\n")   

def getOverlapRegion(l_region, r_region):
    (l_chr,l_start, l_end, l_direct) = (l_region.chr,l_region.start, l_region.end, l_region.direction)
    (r_chr, r_start, r_end, r_direct) = (r_region.chr,r_region.start, r_region.end, r_region.direction)
    l_reads = l_region.for_reads+l_region.rev_reads
    r_reads = r_region.for_reads+r_region.rev_reads
    read_ratio = float(l_reads)/(r_reads)
    if l_chr !=r_chr or l_direct==r_direct or read_ratio<0.1 or read_ratio>10:
        return ""
    else:
        start = max(l_start, r_start)
        end = min(l_end, r_end)
        if start<end:
            return l_chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(end-start)+'\t'+str(l_reads+r_reads)+"\t"+l_chr+":"+str(start)+"-"+str(end)
        else:
            return ""

def addRegionToDic(region_dic, s_id, region):
    chr, start, end = (region.chr, region.start, region.end)
    if chr not in region_dic:
        region_dic[chr]=[]
    target_region = None
    # Search whether the region is already in the list
    for s in region_dic[chr]:
        if getOverlapRegion(s["reg"],region)!="":
            s["ids"].append(s_id)
            target_region = s["reg"]
            break
    if target_region==None:
        region_dic[chr].append({"ids":[s_id],"reg":region})
    return target_region

def siteDic2RegionDic(site_dic):
    region_dic={}
    for i,s in site_dic.items():
        (l_chr, l_loc, l_direct) = (s.l_chr, s.l_loc, s.l_direct)
        l_id=l_chr+"_"+str(l_loc)+"_"+l_direct
        (r_chr, r_loc, r_direct) = (s.r_chr, s.r_loc, s.r_direct)
        r_id=r_chr+"_"+str(r_loc)+"_"+r_direct
        if l_chr not in region_dic:
            region_dic[l_chr]={}
        if r_chr not in region_dic:
            region_dic[r_chr]={}
        region_dic[l_chr][l_id]=(l_chr, l_loc, l_direct, s.for_read_num, s.rev_read_num)
        region_dic[r_chr][r_id]=(r_chr, r_loc, r_direct, s.for_read_num, s.rev_read_num)
    return region_dic

def getGenomeSizeDic(genomeSize_file):
    genomeSizeDic={}
    with open(genomeSize_file, 'r') as SIZE:
        for l in SIZE:
            contents = l.strip().split("\t")
            chr = contents[0]
            size = int(contents[1])
            genomeSizeDic[chr]=size
    return genomeSizeDic

def annoSiteByBigWig2(site_dic, genomeSize_file, bw_file, outfile, cutoff=30):
    bw = pyBigWig.open(bw_file)
    site_region_dic = siteDic2RegionDic(site_dic)
    genomeSize_dic= getGenomeSizeDic(genomeSize_file)
    l_region_dic={}
    r_region_dic={}
    filtered_site_dic={}
    interaction_dic = {}
    for chr,region_dic in site_region_dic.items():
        if chr in genomeSize_dic:
            print(chr)
            size = genomeSize_dic[chr]
            temp_ls=bw.values(chr, 0, size-1)
            for reg_id, values in region_dic.items():
                pos = int(values[1])
                direct= values[2]
                region_dic[reg_id]=findRegion2(temp_ls,chr, pos, direct, cutoff)
            temp_ls=None
    with open(outfile, 'w') as OUTPUT:
        for i,s in site_dic.items():
            s_id = s.getID()
            [l_reg_id, r_reg_id] = s.getID().split("#")
            l_chr = s.l_chr
            r_chr = s.r_chr
            l_reg = site_region_dic[l_chr][l_reg_id]
            r_reg = site_region_dic[r_chr][r_reg_id]
            overlap_region=""
            if l_chr in genomeSize_dic and r_chr in genomeSize_dic:
                if l_reg.isValidated() and r_reg.isValidated():
                    overlap_region = getOverlapRegion(l_reg,r_reg)
                    if overlap_region=="":
                        #addRegionToDic(l_region_dic, s_id, l_reg)
                        #addRegionToDic(r_region_dic, s_id, r_reg)
                        filtered_site_dic[s_id] = {"l":l_reg, "r":r_reg}
                    else:
                        # Self connect
                        contents = overlap_region.split("\t")
                        chr=contents[0]
                        start = contents[1]
                        end = contents[2]
                        length = contents[3]
                        reads= contents[4]
                        OUTPUT.write(chr+"\t"+start+"\t"+end+"\t"+s_id+"\t"+reads+"\t"+l_reg.getID()+"\t"+str(l_reg.for_read_num)+"/"+str(l_reg.rev_read_num)+"\t"+r_reg.getID()+"\t"+str(r_reg.for_read_num)+"/"+str(r_reg.rev_read_num)+"\t"+overlap_region+"\n")
    # Connect the multiple /singlon eccDNAs
    interact_pair_ls = []
    for s_id, s_obj in filtered_site_dic.items():
        s_r_reg = s_obj["r"]  # upstream
        s_l_reg = s_obj["l"]  # upstream
        for ss_id, ss_obj in filtered_site_dic.items():
            if s_id != ss_id:
                ss_l_reg = ss_obj["l"]
                ss_r_reg = ss_obj["r"]
                
                overlap = getOverlapRegion(s_r_reg, ss_l_reg)
                if (overlap!="" and abs(int(overlap.split("\t")[3]))>2):
                    # s_id -> ss_id
                    pair_str = s_id+"\t"+ss_id
                    if pair_str not in interact_pair_ls:
                        interact_pair_ls.append(pair_str)
                
                overlap = getOverlapRegion(ss_r_reg, s_l_reg)
                if (overlap!="" and abs(int(overlap.split("\t")[3]))>2):
                    # ss_id -> s_id
                    pair_str = ss_id+"\t"+s_id
                    if pair_str not in interact_pair_ls:
                        interact_pair_ls.append(pair_str)
    
    # for i in filtered_site_dic:
    #     cs = filtered_site_dic[i]
    #     s_r_reg = cs["r"]
    #     s_r_chr = s_r_reg.chr
    #     if s_r_chr in l_region_dic:
    #         l_reg_ls = l_region_dic[s_r_chr]
    #         for l_reg in l_reg_ls:
    #             # site->region
    #             overlap = getOverlapRegion(l_reg["reg"], s_r_reg)
    #             if (overlap!="" and abs(int(overlap.split("\t")[3]))>2):
    #                 # macthed
    #                 if i not in interaction_dic:
    #                     interaction_dic[i]=[]
    #                 interaction_dic[i]=list(set(interaction_dic[i]+l_reg["ids"]))
    #     s_l_reg = cs["l"]
    #     s_l_chr = s_l_reg.chr
    #     for s_l_chr in r_region_dic:
    #         r_reg_ls = r_region_dic[s_l_chr]
    #         for r_reg in r_reg_ls:
    #             # region->site
    #             overlap = getOverlapRegion(r_reg['reg'], s_l_reg)
    #             if (overlap!="" and abs(int(overlap.split("\t")[3]))>2):
    #                 # macthed
    #                 for s_id in r_reg["ids"]:
    #                     if (s_id not in interaction_dic):
    #                         interaction_dic[s_id] = []
    #                     interaction_dic[s_id]=list(set(interaction_dic[s_id]+[i]))
    # return interaction_dic
    return interact_pair_ls

def annoSiteByBigWig3(site_dic, genomeSize_file, bw_file, outfile, bed_file, cutoff=30):
    bw = pyBigWig.open(bw_file)
    site_region_dic = siteDic2RegionDic(site_dic)
    genomeSize_dic= getGenomeSizeDic(genomeSize_file)
    l_region_dic={}
    r_region_dic={}
    filtered_site_ls=[]
    for chr,region_dic in site_region_dic.items():
        if chr in genomeSize_dic:
            print(chr)
            size = genomeSize_dic[chr]
            temp_ls=bw.values(chr, 0, size-1)
            for reg_id, values in region_dic.items():
                pos = int(values[1])
                direct= values[2]
                for_reads=values[3]
                rev_reads=values[4]
                region_dic[reg_id]=findRegion2(temp_ls,chr, pos, direct, for_reads, rev_reads, cutoff)
            temp_ls=None
    with open(outfile, 'w') as OUTPUT, open(bed_file,"w") as BED:
        for i,s in site_dic.items():
            s_id = s.getID()
            [l_reg_id, r_reg_id] = s.getID().split("#")
            l_chr = s.l_chr
            r_chr = s.r_chr
            l_reg = site_region_dic[l_chr][l_reg_id]
            r_reg = site_region_dic[r_chr][r_reg_id]
            overlap_region=""
            if l_chr in genomeSize_dic and r_chr in genomeSize_dic:
                if l_reg.isValidated() and r_reg.isValidated():
                    # Test if it is self-connection
                    overlap_region = getOverlapRegion(l_reg,r_reg)
                    if overlap_region=="":
                        #addRegionToDic(l_region_dic, s_id, l_reg)
                        #addRegionToDic(r_region_dic, s_id, r_reg)
                        filtered_site_ls.append({"id":s_id, "l":l_reg, "r":r_reg})
                    else:
                        contents = overlap_region.split("\t")
                        chr=contents[0]
                        start = contents[1]
                        end = contents[2]
                        length = contents[3]
                        reads= contents[4]
                        BED.write(chr+"\t"+start+"\t"+end+"\t"+s_id+"\n")
                        OUTPUT.write(s_id+"\t"+l_reg.getID()+"\t"+str(l_reg.coverage)+"\t"+r_reg.getID()+"\t"+str(r_reg.coverage)+"\t"+overlap_region+"\n")
    # Connect the multiple /singlon eccDNAs
    interact_pair_ls = []
    for i in range(len(filtered_site_ls)-1):
        s_obj = filtered_site_ls[i]
        s_id =  s_obj["id"]
        s_r_reg = s_obj["r"]  # upstream
        s_l_reg = s_obj["l"]  # upstream
        for k in range(i+1, len(filtered_site_ls)):
            ss_obj = filtered_site_ls[k]
            ss_id = ss_obj["id"]
            ss_l_reg = ss_obj["l"]
            ss_r_reg = ss_obj["r"]
            
            overlap = getOverlapRegion(s_r_reg, ss_l_reg)
            if (overlap!="" and abs(int(overlap.split("\t")[3]))>2):
                # s_id -> ss_id
                pair_str = s_id+"\t"+ss_id
                if pair_str not in interact_pair_ls:
                    interact_pair_ls.append(pair_str)
                
            overlap = getOverlapRegion(ss_r_reg, s_l_reg)
            if (overlap!="" and abs(int(overlap.split("\t")[3]))>2):
                # s_id -> ss_id
                pair_str = ss_id+"\t"+s_id
                if pair_str not in interact_pair_ls:
                    interact_pair_ls.append(pair_str)
    return interact_pair_ls

def printInteractionList(interact_pair_ls, network_file, node_file, ecc_file):
    nodes_dic = {}
    index = 0 # The number of interaction pair
    with open(network_file, 'w') as NETWORK:
        for pair in interact_pair_ls:
            l_s_id = pair.split("\t")[0]
            r_s_id = pair.split("\t")[1]
            if l_s_id not in nodes_dic:
                nodes_dic[l_s_id] = index
                index+=1
            if r_s_id not in nodes_dic:
                nodes_dic[r_s_id] = index
                index+=1
            NETWORK.write(pair+"\n")
    # Build the graph and find the eccDNAs
    g = Graph(index) 
    for pair in interact_pair_ls:
        l_s_id = pair.split("\t")[0]
        l_s_index = nodes_dic[l_s_id]
        r_s_id = pair.split("\t")[1]
        r_s_index = nodes_dic[r_s_id]
        g.addEdge(l_s_index, r_s_index)

    # print node file and networks files
    with open(node_file, 'w') as NODE:
        for id, index in nodes_dic.items():
            NODE.write(str(index)+"\t"+str(id)+"\n")
    g.printSCCs(ecc_file)

def printInteractionDic(inter_dic, network_file, node_file, ecc_file):
    nodes_dic = {}
    index = 0
    with open(network_file, 'w') as NETWORK:
        for id in inter_dic:
            if id not in nodes_dic:
                nodes_dic[id]=index
                index+=1
            target_ids = inter_dic[id]
            for t_id in target_ids:
                NETWORK.write(id+"\t"+t_id+"\n")
                if t_id not in nodes_dic:
                    nodes_dic[t_id]=index
                    index+=1
    # Build the graph and find the eccDNAs
    g = Graph(index) 
    for id in inter_dic:
        target_ids = inter_dic[id]
        s_index =  nodes_dic[id]
        for t_id in target_ids:
            e_index = nodes_dic[t_id]
            g.addEdge(s_index, e_index)
    # print node file and networks files
    with open(node_file, 'w') as NODE:
        for id, index in nodes_dic.items():
            NODE.write(str(index)+"\t"+str(id)+"\n")
    g.printSCCs(ecc_file)

def reconstruct_multi_element_ecc(network_file, mut_ecc_file):
    node_dic=[]
    with open(network_file,'r') as NODE:
        for line in NODE:
            contents = line.strip().split("\t")
            if len(contents)==2:
                node_dic[contents[0]]= contents[1]
    with open(ecc_file,"r") as ECC:
        for line in ECC:
            contents = line.strip().strip(',').split(",")
            if len(contents)>1:
                pass

# chr1_201748285_+#chr1_201749248_+,chr1_201750949_-#chr1_201748782_-
# chr13_82392342_-#chr13_82390110_-,chr13_82388097_+#chr13_82386482_+
def getRegionStr(pair_ls):
    c = pair_ls
    a = c[0].split("#")[1].split("_")
    b = c[1].split("#")[0].split("_")
    return a[0]+":"+a[1]+"-"+b[1]
        
    

def findeccDNA(net_file, node_file, ecc_file, outfile):
    net = pd.read_csv(net_file, sep="\t",header=None)
    net.columns=['b_s_id', 'e_s_id']
    node = pd.read_csv(node_file, sep="\t",header=None)
    node.columns = ['index','s_id']
    eccDNA_ls= []
    with open(ecc_file, 'r') as ECC:
        for line in ECC:
            # if line.strip()=="1663,1409,":
            #     print(line)
            contents =list(map(int, line.strip().strip(',').split(",")))
            if len(contents)>1:
                sel_nodes= node[node['index'].isin(contents)]["s_id"].to_list()
                if len(sel_nodes)>1:
                    sel_net = net[net['b_s_id'].isin(sel_nodes) & net['e_s_id'].isin(sel_nodes)]
                    pair_ls = list(zip(sel_net['b_s_id'].tolist(), sel_net['e_s_id'].tolist()))
                    for i in range(0, len(pair_ls)-1):
                        eccDNA=[pair_ls[i]]
                        for k in range(i+1, len(pair_ls)):
                            if (eccDNA[-1][1]==pair_ls[k][0]):
                                eccDNA.append(pair_ls[k])
                                if (pair_ls[k][1]==eccDNA[0][0]):
                                    print(len(eccDNA))
                                    eccDNA_ls.append(list(eccDNA))
                                    break
                            elif (eccDNA[0][0]==pair_ls[k][1]):
                                eccDNA=[pair_ls[k]]+eccDNA
                                if (pair_ls[k][1]==eccDNA[0][0]):
                                    print(len(eccDNA))
                                    eccDNA_ls.append(list(eccDNA))
                                    break
                        
    with open(outfile, 'w') as OUTPUT:
        for i in eccDNA_ls:
            line_ls =[]
            for k in i:
                line_ls.append(getRegionStr(k))
            OUTPUT.write("\t".join(line_ls)+"\n")


# junc_site_dic=getJunSiteDic(sample_bam)
# cons_fasta =out_dir+"/cons.fa"
# JunSiteDicToFasta(junc_site_dic, cons_fasta)
# blast_out = BLASTNConsenseSeq("blastn","makeblastdb","test",ref_fa,cons_fasta,out_dir)

# #blast_out="/media/hqyone/2tb/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704/test_tRNA_blast_out.tab"
# (site_hit_dic, Site2_dic) = ParseBlastResult(blast_out,junc_site_dic)
# #bed_file="/home/hqyone/mnt/2t/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak"
# bigwig_file=test_dir+ "/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704/macs/sample1.bigwig"
# genome_size_file = "/media/hqyone/2tb/eccDNA/genome/chrom_size/hg38.chrom.sizes"
# #refined_bed_file="/home/hqyone/mnt/2t/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.refine"

# #refineBedFile(bed_file, bigwig_file, refined_bed_file)
# out_file=test_dir+"/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.output"
# bed_file=test_dir+"/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.output.bed"
# net_file=test_dir+"/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.network"
# node_file =test_dir+"/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.node"
# ecc_file =test_dir+"/eccDNA/data-raw/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.ecc"

# #int_dic = annoSiteByBigWig(Site2_dic, bigwig_file, out_file)
# #int_dic = annoSiteByBigWig2(Site2_dic, genome_size_file,bigwig_file, out_file)
# #printInteractionDic(int_dic,net_file, node_file, ecc_file)


# int_ls= annoSiteByBigWig3(Site2_dic, genome_size_file,bigwig_file, out_file, bed_file)
# printInteractionList(int_ls,net_file, node_file, ecc_file)

#annoSiteByBed(Site2_dic, refined_bed_file, out_file)
#print(site_hit_dic)











