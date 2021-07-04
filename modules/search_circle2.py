'''
Author: Quanyuan(Leo) He
Email: hqyone@gmail.com
Insititute: Hunan Normal Univeristy
Date: 1969-12-31 17:00:00
LastEditTime: 2021-04-20 20:59:23
LastEditors: Quanyuan(Leo) He
Description: 
FilePath: /eccDNA/code/source/search_circle2.py
License: The MIT License (MIT)
'''
import pandas as pd

def isCrossChrSite(s_id):
    return s_id.split("#")[0].split("_")[0] != s_id.split("#")[1].split("_")[0]

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
                
  
# Thanks to Divyanshu Mehta for contributing this code                
            
            
test_dir = "/media/hqyone/2tb1"
net_file=test_dir+"/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.network"
#net_file=test_dir+"/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/test.network"
node_file =test_dir+"/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.node"
ecc_file =test_dir+"/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.ecc"
ecc_out_file =test_dir+"/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.eccout.csv"

findeccDNA(net_file, node_file, ecc_file, ecc_out_file)