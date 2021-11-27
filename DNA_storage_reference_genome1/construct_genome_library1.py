# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:16:50 2020

@author: Jitao Zhang
"""
import struct
import sys
sys.path.append("../../../")
import LZW.compress_short
# import DNAstorage_reference_lookup.genome
import genome
import numpy as np
import pandas as pd
import time
import re
from openpyxl import load_workbook, Workbook
import math
import os


def decoded_method(mode=None, genome_nt3_dict=None, seed=None):
    decode_dict = {}

    decode_dict_old = {'AAA': 31, 'TTT': 32, 'AAT': 30, 'ATT': 33,
                       'ATA': 29, 'TAT': 34, 'GAA': 28, 'TTC': 35,
                       'CAA': 27, 'TTG': 36, 'TAA': 26, 'TTA': 37,
                       'AAG': 25, 'CTT': 38, 'AGA': 24, 'TCA': 39,
                       'TGA': 23, 'TCT': 40, 'CAT': 22, 'ATG': 41,
                       'AAC': 21, 'GTT': 42, 'ATC': 20, 'GAT': 43,
                       'ACA': 19, 'TGT': 44, 'AGT': 18, 'ACT': 45,
                       'CCA': 17, 'TGG': 46, 'GTA': 16, 'TAC': 47,
                       'TAG': 15, 'CTA': 48, 'GGA': 14, 'TCC': 49,
                       'CAG': 13, 'CTG': 50, 'GCA': 12, 'TGC': 51,
                       'ACC': 11, 'AGC': 52, 'GCT': 10, 'AGG': 53,
                       'GGT': 9, 'CCT': 54, 'GAG': 8, 'CTC': 55,
                       'CAC': 7, 'GTG': 56, 'GAC': 6, 'GTC': 57,
                       'CGA': 5, 'TCG': 58, 'ACG': 4, 'CGT': 59,
                       'GCC': 3, 'GGC': 60, 'CCC': 2, 'GGG': 61,
                       'CCG': 1, 'CGG': 62, 'CGC': 0, 'GCG': 63}


    if mode == 'sublib':
        if genome_nt3_dict.keys() == decode_dict_old.keys():
            print("genome_nt3_dict's keys is equal to decode_dict_old's keys()")
            return decode_dict_old
        else:
            print("genome_nt3_dict's keys is not equal to decode_dict_old's keys()")
            key_list = list(genome_nt3_dict.keys())
            for index, value in enumerate(decode_dict_old.values()):
                decode_dict[key_list[index]] = value

    if mode == 'random':
        np.random.seed(seed=seed)
        random_list = np.random.choice(range(0, 64), 64, 0)
        for index, key in enumerate(decode_dict_old.keys()):
            decode_dict[key] = random_list[index]
    elif mode == 'reversed':
        random_list = [i for i in decode_dict_old.values()]
        random_list = random_list[::-1]
        for index, key in enumerate(decode_dict_old.keys()):
            decode_dict[key] = random_list[index]
    else:
        decode_dict = decode_dict_old

    # print("decode_dict is:",decode_dict)
    return decode_dict


def calculate_gc_content(self, encoded_data_list):
    encoded_data_seq = {'Y_channel': '', 'Cb_channel': '', 'Cr_channel': ''}
    two_nts_dict = {'AA': 0, 'AT': 0, 'AC': 0, 'AG': 0, 'TA': 0, 'TT': 0, 'TC': 0, 'TG': 0, 'CA': 0, 'CT': 0, 'CC': 0,
                    'CG': 0, 'GA': 0, 'GT': 0, 'GC': 0, 'GG': 0}
    nt_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    for index_of_channel, channel in enumerate(encoded_data_list):

        for data_matrix in encoded_data_list[channel]:
            seq_of_data_matrix = data_matrix[0]
            for difference in data_matrix[1]:
                seq_of_data_matrix += self.encode_method(difference)
            encoded_data_seq[channel] += seq_of_data_matrix

    for channel in encoded_data_seq:
        for seq in encoded_data_seq[channel]:
            for nt in seq:
                nt_dict[nt] += 1

        for index_of_two_nts in range(len(encoded_data_seq[channel]) - 1):
            two_nts_dict[encoded_data_seq[channel][index_of_two_nts:index_of_two_nts + 2]] += 1

    nt_ratio = [nt / sum(nt_dict.values()) for nt in nt_dict.values()]

    two_nts_ratio = [nt2 / sum(two_nts_dict.values()) for nt2 in two_nts_dict.values()]

    print("nt_ratio is: ", nt_ratio)
    print("the sum of nt is:", sum(nt_dict.values()))
    print("two_nts_ratio is: ", two_nts_ratio)
    print("the sum of 2nt is:", sum(two_nts_dict.values()))
    return None


def construct_dict(path, kmer, size, seed):
    gn = genome.genome(path)
    if size:
        lib = gn.construct_library(kmer, [size, seed])

    else:
        lib = gn.construct_library(kmer)
    # decoded_dict = gn._get_nt_ratio(lib,1,1)[2]
    print("power is {}".format(power), "library{} complete".format(kmer))

    print("the lenght of gn.library is :" + " " + str(len(lib)))

    gn_dict = {}
    count = 0
    cycle_init = time.time()

    for segment in lib:
        if segment[:4] not in gn_dict:
            gn_dict[segment[:4]] = []
        gn_dict[segment[:4]].append(segment[4:])

        count += 1
        if count % 10000 == 0:
            print("\r{:.2%} has completed,this cycle cost: {:.2f}s".format(count / len(lib), time.time() - cycle_init),
                  end='')
            cycle_init = time.time()

    # for key,value in gn_dict.items():
    #    gn_dict[key] = np.array(value)

    genome_nt_dict, genome_nt3_dict, genome_nt4_dict = gn.get_genome_dict_ratio(gn_dict)
    genome_data = [gn_dict, genome_nt_dict, genome_nt3_dict, genome_nt4_dict]
    print("gn.lib completed")
    return genome_data


if __name__ == '__main__':
    genome_range = [15, 20]
    genome_list = pd.read_excel("E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\kl-test\genome_list.xlsx",
                                sheet_name="genome-size").iloc[genome_range[0]:genome_range[1], 0].values
    seed = 2020
    kmer_list = [16, 52]
    power_list = [7, 8, 9, 10, 11]
    
    
    for genome_name in genome_list:
        path = r"./lena512_16_16_matrix/kl-test/download-genome/multi-genome/{}/{}".format(genome_name + ".fna",
                                                                                           genome_name + ".fna")

        # writer = pd.ExcelWriter(r'E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\color-bmp\result\genome-result\nt4_nt3_info.xlsx',engine='openpyxl')

        for power in power_list:
            size = 4 ** power
            mode = 'random'
            for kmer in kmer_list:
                genome_data = construct_dict(path, kmer, size, seed)
                if not os.path.exists("G:\work_in_cngb\genome-size\genome_dict\{}".format(genome_name)):
                    os.mkdir("G:\work_in_cngb\genome-size\genome_dict\{}".format(genome_name))
                np.save(r"G:\work_in_cngb\genome-size\genome_dict\{}\genome-{}-seed_{}-power_{}-bioseq.npy".format(
                    genome_name, kmer, seed, power), genome_data)
    