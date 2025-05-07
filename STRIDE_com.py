#-*- coding:utf-8 -*-
import sys
import numpy as np
import os
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
import re



def check_seq(seq1):
    n1 = seq1.count("A")
    n2 = seq1.count("T") + seq1.count("U")
    n3 = seq1.count("C")
    n4 = seq1.count("G")
    if len(seq1) == n1 + n2 + n3 + n4:
        return(1);
    else:
        return(0);


def read_fa(fa_file, dict1):
    """
    >ENST00000622302.2|ENSG00000277507.2|KB-226F1.4|unprocessed_pseudogene
    GAGCTGTTTCCGTTCCTCTGCCCGCCATGCCGTTCCTAGAGCTGCACACGAATTTCCCCGCCAACCGAGTGCCCGCGG
    """
    dict2 = {}
    f = open(fa_file, 'r')
    i = 0
    r = 0
    sid, name, sequ = "", "NULL", ""
    for line in f:
        i = i + 1
        line = line.strip('\n')
        if i % 2 == 1:
            r = r + 1
            sent1 = line[1:].split(".")
            sid = sent1[0]
            name = "NULL"
            sequ = ""
        else:
            sequ = line
            sequ = sequ.upper().replace("U", "T")
            if sid in dict1:
                name = dict1[sid]
            dict2["query_" + str(r)] = [sequ, len(sequ), sid, name]
    f.close()
    return dict2


def read_trx_info(trx_file):
    """
    >ENST00000622302.2|ENSG00000277507.2|KB-226F1.4|unprocessed_pseudogene
    GAGCTGTTTCCGTTCCTCTGCCCGCCATGCCGTTCCTAGAGCTGCACACGAATTTCCCCGCCAACCGAGTGCCCGCGG
    """
    dict1, dict2 = {}, {}
    f = open(trx_file, 'r')
    i = 0
    #sid, name, sequ = "", "", ""
    for line in f:
        i = i + 1
        line = line.strip('\n')
        if line[0] == ">":
            sent1 = line[1:].split("|")
            sent2 = sent1[0].split(".")
            sent3 = sent1[1].split(".")
            #sent4 = sent3[1].split(".")
            trx_id, gene_name, gene_id = sent2[0], sent1[2], sent3[0]
            dict1[trx_id] = gene_id
            if gene_id in dict2:
                dict2[gene_id].append(trx_id)
            else:
                dict2[gene_id] = [trx_id]
    f.close()
    return dict1, dict2


def read_trx_info2(trx_file):
    """
    ENST00000000233.9       ARF5=ENSG00000004059.10 protein_coding  1103    1-221,222-302,303-412,413-484,485-610,611-1103  155     695
    ENST00000000412.7       M6PR=ENSG00000003056.7  protein_coding  2756    1-468,469-645,646-812,813-922,923-1053,1054-1180,1181-2756    470      1301
    """
    dict1, dict2 = {}, {}
    f = open(trx_file, 'r')
    i = 0
    #sid, name, sequ = "", "", ""
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent1 = line.split("\t")
        sent2 = sent1[0].split(".")
        sent3 = sent1[1].split("=")
        sent4 = sent3[1].split(".")
        trx_id, gene_name, gene_id = sent2[0], sent3[0], sent4[0]
        dict1[trx_id] = gene_id
        if gene_id in dict2:
            dict2[gene_id].append(trx_id)
        else:
            dict2[gene_id] = [trx_id]
    f.close()
    return dict1, dict2


def read_kmer(kmer_file):
    """
    TGTGTGTGTGTGTGT    30298    0.000142963518511025
    """
    dict1, dict2 = {}, {}
    f = open(kmer_file, 'r')
    i = 0
    #sid, name, sequ = "", "", ""
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent1 = line.split("\t")
        dict1[sent1[0]] = [int(sent1[1]), float(sent1[2])]
    f.close()
    return dict1


def read_shape(shape_file):
    """
    shape_file       -- A file of shape data
    ENST00000341426 3235    7.714   NULL    NULL    NULL    NULL    NULL
    """
    dict1, dict2 = {}, {}
    f = open(shape_file, 'r')
    i = 0
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent = line.split('\t')
        sent1 = sent[0].split('.')
        #if cal_coverage(sent[3:]) >= 0.5:
        dict1[sent1[0]] = sent[3:]
        #dict2[sent1[0]] = [int(sent[1]), sent[2]]
    #plen = i
    f.close()
    return dict1


def read_shape2(shape_file1, shape_file2):
    """
    shape_file       -- A file of shape data
    ENST00000341426 3235    7.714   NULL    NULL    NULL    NULL    NULL
    """
    dict1, dict2 = {}, {}
    f = open(shape_file1, 'r')
    i = 0
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent = line.split('\t')
        sent1 = sent[0].split('.')
        #if cal_coverage(sent[3:]) >= 0.5:
        dict1[sent1[0]] = sent[3:]
        #dict2[sent1[0]] = [int(sent[1]), sent[2]]
    #plen = i
    f.close()
    f = open(shape_file2, 'r')
    i = 0
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent = line.split('\t')
        sent1 = sent[0].split('.')
        if sent1[0] in dict1:
            shape1 = dict1[sent1[0]]
            shape2 = sent[3:]
            shape_list = []
            for j in range(len(shape1)):
                if shape1[j] != "NULL" and shape2[j] != "NULL":
                    shape_list.append(str(float(shape1[j]) - float(shape2[j])))
                else:
                    shape_list.append("NULL")
            dict2[sent1[0]] = shape_list
        #if cal_coverage(sent[3:]) >= 0.5:
        #dict1[sent1[0]] = sent[3:]
        #dict2[sent1[0]] = [int(sent[1]), sent[2]]
    #plen = i
    f.close()
    return dict2


def filter_Tm(sdict, length1, kmer_dict):
    #sdict key : sequence,length,trx_id,gene_id
    pdict = {}
    for ke in sdict:
        sequ = sdict[ke][0]
        pdict[ke] = {}
        for i in range(len(sequ) - length1 + 1):
            seq1 = sequ[i:(i+length1)]
            flag = 0
            list1 = ["AAAAA", "CCCCC", "GGGGG", "TTTTT"]
            for ke2 in list1:
                if ke2 in seq1:
                    flag = 1
                    break
            if flag == 1:
                continue
            for j in range(len(seq1) - 15 + 1):
                seq2 = seq1[j:(j+15)]
                if seq2 in kmer_dict:
                    flag = 1
                    break
            if flag == 1:
                continue
            n1 = seq1.count("C") + seq1.count("G")
            n2 = len(seq1)
            gc_content = 1.0 * n1 / n2
            if gc_content >= 0.25 and gc_content <= 0.75:
                pdict[ke][i] = [gc_content, abs(gc_content - 0.5), seq1]
    return pdict


def filter_spec(sdict, pdict, path1, trx_file, cut_off):
    #sdict key -> sequence,length,trx_id,gene_id
    #pdict key -> position -> gc_content, abs(gc_content - ave), probe_region
    pdict2 = pdict
    probe_file, blast_out = path1 + "tmp1_prespec.fa", path1 + "tmp1_blast.txt"
    fw1 = open(probe_file, 'w')
    #fw2 = open(tmpfile2, 'w')
    for sid in sdict:
        sequence, length, trx_id, gene_id = sdict[sid][0], sdict[sid][1], sdict[sid][2], sdict[sid][3]
        for ke in pdict[sid]:
            outstr1 = ">" + sid + "|" + str(ke) + "|" + gene_id + "|" + trx_id + "\n" + pdict[sid][ke][2] + "\n"
            fw1.write(outstr1)
    fw1.close()
    cmd = "blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 " + trx_file + " " + probe_file + " " + blast_out
    os.system(cmd)
    f1 = open(blast_out, 'r')
    for line in f1:
        line = line.strip('\n')
        sent1 = line.split("\t")
        if sent1[0].isnumeric() and int(sent1[0]) >= cut_off and sent1[8] == "+":
            sent2 = sent1[9].split("|")
            sid, ke, gene_id = sent2[0], int(sent2[1]), sent2[2]
            sent3 = sent1[13].split("|")
            sent4 = sent3[1].split(".")
            gene_id2 = sent4[0]
            if gene_id != gene_id2:
                #import pdb; pdb.set_trace()
                if ke in pdict2[sid]:
                    del pdict2[sid][ke]
                #else:
                    #print(sid + "\t" + str(ke) + "\n")
    f1.close()
    return pdict2
    

def filter_hybrid(sdict, pdict, path1, cut_off, length1):
    #sdict key -> sequence,length,trx_id,gene_id
    #pdict key -> position -> gc_content, abs(gc_content - ave), probe_region
    pdict2 = pdict
    hdict = {}
    probe_file, rev_probe_file, blast_out = path1 + "tmp1_reg1.fa", path1 + "tmp1_reg2.fa", path1 + "tmp1_blast2.txt"
    fw1 = open(probe_file, 'w')
    fw2 = open(rev_probe_file, 'w')
    for sid in pdict:
        for ke in pdict[sid]:
            seq1 = pdict[sid][ke][2]
            seq2 = seq1[::-1].replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c")
            seq2 = seq2.upper()
            outstr1 = ">" + sid + "|" + str(ke) + "\n" + seq1 + "\n"
            fw1.write(outstr1)
            outstr2 = ">" + sid + "|" + str(ke) + "\n" + seq2 + "\n"
            fw2.write(outstr2)
    fw1.close()
    fw2.close()
    cmd = "blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 " + rev_probe_file + " " + probe_file + " " + blast_out
    os.system(cmd)
    f1 = open(blast_out, 'r')
    for line in f1:
        line = line.strip('\n')
        sent1 = line.split("\t")
        if sent1[0].isnumeric() and int(sent1[0]) >= cut_off and sent1[8] == "+":
            sid1, sid2 = sent1[9], sent1[13]
            if sid1 == sid2:
                sent2 = sent1[9].split("|")
                sid, ke = sent2[0], int(sent2[1])
                del pdict2[sid][ke]
            else:
                if sid1 in hdict:
                    hdict[sid1].append(sid2)
                else:
                    hdict[sid1] = [sid2]
    f1.close()
    for ke in hdict:
        sent2 = ke.split("|")
        sid, ke = sent2[0], int(sent2[1])
        if length(hdict[ke]) >= 3:
            del pdict2[sid][ke]
    odict = probe_overlap(sdict, pdict2, length1)
    dudict = {}
    for ke in odict:
        if ke in dudict:
            sent2 = ke.split("|")
            sid, ke = sent2[0], int(sent2[1])
            del pdict2[sid][ke]
        if ke in hdict:
            sent1 = hdict[ke]
            for ke2 in sent1:
                dudict[ke2] = 1
    return pdict2, hdict


def probe_overlap(sdict, pdict, length1):
    odict = {}
    for sid in pdict:
        list1 = [k for k in pdict[sid].keys()]
        list1.sort()
        length2 = sdict[sid][1]
        list2 = [0 for i in range(length2)]
        #import pdb; pdb.set_trace()
        for k in list1:
            for j in range(length1):
                list2[k+j] = list2[k+j] + 1
        for ke in pdict[sid]:
            sid2 = sid + "|" + str(ke)
            n1 = max(list2[ke:(ke+length1)])
            odict[sid2] = n1
    dict2 = {k: v for k, v in sorted(odict.items(), key=lambda item: item[1])}
    return dict2


def cal_mean(sent1, ke, length1):
    n1, n2 = 0.0, 0.0
    for i in range(ke, ke+length1):
        if sent1[i] != "NULL":
            n1 = n1 + float(sent1[i])
            n2 = n2 + 1
    n3, n4 = 0, 0
    for i in range(len(sent1)):
        if sent1[i] != "NULL":
            n3 = n3 + float(sent1[i])
            n4 = n4 + 1
    if n2 >= length1 / 4:
        return n1 / n2
    else:
        return n3 / n4


def cal_probe_shape(sdict, pdict, shapedict, length1):
    for sid in pdict:
        trx_id, gene_id = sdict[sid][2], sdict[sid][3]
        sent1 = shapedict[trx_id]
        for ke in pdict[sid]:
            pdict[sid][ke][1] = cal_mean(sent1, ke, length1)
    return pdict


def define_state(sdict, pdict, rlength, length1):
    pdict2 = {}
    for sid in sdict:
        sequence, length, trx_id, gene_id = sdict[sid][0], sdict[sid][1], sdict[sid][2], sdict[sid][3]
        if sid in pdict:
            pdict2[sid] = {}
            for k in pdict[sid]:
                r = int(k / rlength)
                if r in pdict2[sid]:
                    pdict2[sid][r].append(k)
                else:
                    pdict2[sid][r] = [k]
    pdict3 = {}
    for sid in pdict2:
        pdict3[sid] = {}
        dict1 = pdict2[sid]
        dict2 = {k: v for k, v in sorted(dict1.items(), key=lambda item: item[0])}
        start_p = 0
        for kk in dict2:
            sent1 = dict2[kk]
            dict3 = {}
            for ke in sent1:
                dict3[ke] = pdict[sid][ke][1]
            dict4 = {k: v for k, v in sorted(dict3.items(), key=lambda item: item[1], reverse = True)}
            #dict4 = {k: v for k, v in sorted(dict3.items(), key=lambda item: item[1])}
            for k in dict4:
                if k >= start_p + 5:
                    start_p = start_p + length1
                    pdict3[sid][k] = pdict[sid][k]
                    break
    return pdict3



def count_probe(pdict):
    num1 = 0
    for sid in pdict:
        for ke in pdict[sid]:
            num1 = num1 + 1
    return num1

def print_probe(sdict, pdict, out_file1, shape_file1, shape_file2, length1):
    num1 = 0
    shapedict1 = read_shape(shape_file1)
    shapedict2 = read_shape(shape_file2)
    fw1 = open(out_file1, 'w')
    for sid in pdict:
        trx_id, gene_id = sdict[sid][2], sdict[sid][3]
        for ke in pdict[sid]:
            #import pdb; pdb.set_trace()
            mean1 = cal_mean(shapedict1[trx_id], ke, length1)
            mean2 = cal_mean(shapedict2[trx_id], ke, length1)
            sent1 = pdict[sid][ke]
            sent2 = [str(ke2) for ke2 in sent1]
            region = sent1[-1]
            probe = region[::-1].replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c")
            probe = probe.upper()
            outstr1 = sid + "\t" + str(ke) + "\t" + trx_id + "\t" + gene_id + "\t" + str(mean1) + "\t" + str(mean2) + "\t" + "\t".join(sent2) + "\t" + probe + "\n"
            fw1.write(outstr1)
    fw1.close()


#/data1/huangwenze/new_z/probe_design/
#python probe_design_v2_com.py HEK293T_hela_dstrx1.fa 60 HEK293T Hela
#python probe_design_v2_com.py hela_HEK293T_dstrx1.fa 60 Hela HEK293T

if __name__ == '__main__':
    seq_file = sys.argv[1]
    length1 = int(sys.argv[2])
    cell_line1 = sys.argv[3]
    cell_line2 = sys.argv[4]

    spec_cutoff = int(length1 / 3 * 2) + 1
    #spec_cutoff = 18
    hybr_cutoff = int(length1 / 3 * 2) + 1
    region_length = 120

    species = "human"
    res_path = "/data1/huangwenze/new_z/probe_design/tmp/"
    data_path = "/data1/huangwenze/new_z/probe_design/data/"
    path1 = res_path + seq_file + "/"
    cmd = "mkdir -p " + path1
    os.system(cmd)
    trx_file = ""
    kmer_file = ""
    if species == "human":
        trx_file = data_path + "hg38_transcriptome_std.fa"
        kmer_file = data_path + "15_mer_frequency_hg38_fa.txt"
    elif species == "mouse":
        trx_file = data_path + "mm10_transcriptome_std.fa"
        kmer_file = data_path + "15_mer_frequency_mm10_fa.txt"

    shape_file1 = data_path + "new_smartSHAPE_0412/" + cell_line1 + "_smartSHAPE.out"
    shape_file2 = data_path + "new_smartSHAPE_0412/" + cell_line2 + "_smartSHAPE.out"

    trx_dict, gene_dict = read_trx_info(trx_file)
    kmer_dict = read_kmer(kmer_file)
    sdict = read_fa(seq_file, trx_dict)
    pdict = filter_Tm(sdict, length1, kmer_dict)
    print(str(count_probe(pdict)) + "\n")
    pdict2 = filter_spec(sdict, pdict, path1, trx_file, spec_cutoff)
    print(str(count_probe(pdict2)) + "\n")
    pdict3, hdict = filter_hybrid(sdict, pdict2, path1, hybr_cutoff, length1)
    print(str(count_probe(pdict3)) + "\n")

    shapedict = read_shape2(shape_file1, shape_file2)
    pdict3 = cal_probe_shape(sdict, pdict, shapedict, length1)

    pdict4 = define_state(sdict, pdict3, region_length, length1)
    print(str(count_probe(pdict4)) + "\n")
    file1 = path1 + seq_file + "_probe1.txt"
    print_probe(sdict, pdict4, file1, shape_file1, shape_file2, length1)
