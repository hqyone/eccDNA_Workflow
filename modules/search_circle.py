'''
Author: Quanyuan(Leo) He
Email: hqyone@gmail.com
Insititute: Hunan Normal Univeristy
Date: 1969-12-31 17:00:00
LastEditTime: 2021-04-16 23:41:25
LastEditors: Quanyuan(Leo) He
Description: 
FilePath: /eccDNA/code/source/search_circle.py
License: The MIT License (MIT)
'''

import os

data=[
    ["A","B"],
    ["B","C"],
    ["C","E"],
    ["G","H"],
    ["H","G"],
    ["C","A"],
    ["A","C"],
    ["M","C"],
    ["N","M"]
]

def removeSelfConnect(data):
    result_data=[]
    selfConnect_ls=[]
    for d in data:
        if len(result_data)==0:
            result_data.append(d)
        else:
            matched=-1
            for i, r in enumerate(result_data):
                if (d[0] == r[1] and d[1]==r[0]):
                    matched=i
                    break
            if matched!=-1:
                selfConnect_ls.append(d)
                del result_data[i]
            else:
                result_data.append(d)
    return (result_data, selfConnect_ls)
                    

def getCircle(data):
    result_data=[]
    for d in data:
        if len(result_data)==0:
            result_data.append(d)
        else:
            matched = 0
            for i,r in enumerate(result_data):
                if d[0]==r[-1]:
                    result_data[i]+=d[1:]
                    matched+=1
                elif d[-1]==r[0]:
                    result_data[i]=d[:-1]+result_data[i]
                    matched+=1
            if matched==0:
                result_data.append(d)
    return result_data


def getAllCircles(data):
    while True:
        result = getCircle(data)
        if result==data:
            return result
        else:
            data=result
            
            
network_file = "/media/hqyone/2tb1/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.network";
network_out_file = "/media/hqyone/2tb1/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.network.out.csv"
network_single_file = "/media/hqyone/2tb1/eccDNA/GHS21010065/Sample_R21009037-1-D21010704/macs2/p_peaks.broadPeak.network.single.csv"
data = []
with open(network_file, 'r') as NETWORK:
    n=0
    for line in NETWORK:
        if n>1000:
            break
        #n+=1
        if line.strip() !="":
            data.append(line.strip().split("\t"))
[data, selfConnect_ls] =  removeSelfConnect(data)
reseult_ls = getAllCircles(data)
with open(network_out_file,'w') as NETWORKOUT:
    for i in reseult_ls:
        NETWORKOUT.write("\t".join(i)+"\n")

with open(network_single_file,'w') as NETWORK_SINGLE_OUT:
    for i in selfConnect_ls:
        NETWORK_SINGLE_OUT.write("\t".join(i)+"\n")

#print("######")
#print(selfConnect_ls)
        
