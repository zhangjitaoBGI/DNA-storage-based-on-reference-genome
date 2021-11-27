# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 11:46:32 2019

@author: 不愿透露姓名的美凌格
"""
import sys
import random
import re
import numpy as np
import pandas as pd
from openpyxl import load_workbook, Workbook
import math
import os
import time

def int_convert_bioseq(number, length=8, low_gc=1):
    bio_seq = ''
    int_to_bin = bin(number)[2:]
    bin_seq = (length - len(int_to_bin)) * '0' + int_to_bin

    for index in range(0, len(bin_seq), 2):
        if bin_seq[index:index + 2] == '00':
            bio_seq += ['C', 'T'][1 - low_gc]
        elif bin_seq[index:index + 2] == '01':
            bio_seq += ['C', 'T'][low_gc]
        elif bin_seq[index:index + 2] == '10':
            bio_seq += ['A', 'G'][1 - low_gc]
        elif bin_seq[index:index + 2] == '11':
            bio_seq += ['A', 'G'][low_gc]
    return bio_seq


def bioseq_convert_int(seq, low_gc=1):
    bin_number = ''
    for element in seq:
        if element == ['C', 'T'][1 - low_gc]:
            bin_number += '00'
        elif element == ['C', 'T'][low_gc]:
            bin_number += '01'
        elif element == ['A', 'G'][1 - low_gc]:
            bin_number += '10'
        elif element == ['A', 'G'][low_gc]:
            bin_number += '11'
    number = int(bin_number, 2)
    return number


class Genome:
    def __init__(self, pth, vitualSize=0):
        '''
        To construct a genome object,
        specify a fasta file contain genome sequences

        Arguments:
            **pth(( --- genome fasta file, set this value None to generate a
            random sequence
            **vitualSize** --- virtual genome size

        Returns:
            *None*
        '''
        print(pth)
        self.pth = pth
        self.data = []
        self.kmer = None
        self.library = {}
        self.sublib = {}
        self.include_N = False
        self.length = 0
        # if self.pth is None:
        #     self.data = ('Virtual Goenome', generateTestSeq.getAseq(vitualSize))
        # else:
        self._loadData()

    def construct_library(self, kmer, sub_lib: list = None):
        '''
        build a library of specify size

        Arguments:
            **self** --- self pointer
            **size** --- segments size

        Return:
            *True* or *False* of this operation
        '''
        if sub_lib:
            size = sub_lib[0]
            seed = sub_lib[1]

        if len(self.data) == 0:
            sys.stderr.write('*** Should load data first ***\n')
            return False
        self.kmer = kmer
        self.library = {}

        offset = 0

        for chromosome in self.data:
            for i in range(len(chromosome[1]) - kmer):
                key = chromosome[1][i:i + kmer]
                if key in self.library:
                    self.library[key].append(i + offset)
                else:
                    self.library[key] = [i + offset, ]
            offset += len(chromosome[1])

        print("The len of genome is:{}".format(offset))

        if sub_lib != None:
            self.sublib = {}
            np.random.seed(seed)
            random_list = np.random.choice(range(len(self.library)), size, 0)
            lib_key = list(self.library.keys())
            random_key = [lib_key[index] for index in random_list]
            for key in random_key:
                self.sublib[key] = self.library[key]
            print("The len of sublib is:{}".format(len(self.sublib)))

            return self.sublib
        else:
            return self.library





    def get_genome_nt_ratio(self, library, output_to_excel=None):
        '''
        get the ratio of each nucleotide or each two nucleotides in the genome
        '''
        length_of_lib = len(library)
        nt_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        nt2_dict = {'AA': 0, 'AT': 0, 'AC': 0, 'AG': 0, 'TA': 0, 'TT': 0, 'TC': 0, 'TG': 0, 'CA': 0, 'CT': 0, 'CC': 0,
                    'CG': 0, 'GA': 0, 'GT': 0, 'GC': 0, 'GG': 0}
        nt3_dict = {'AAA': 0, 'TTT': 0, 'AAT': 0, 'ATT': 0,
                    'ATA': 0, 'TAT': 0, 'GAA': 0, 'TTC': 0,
                    'CAA': 0, 'TTG': 0, 'TAA': 0, 'TTA': 0,
                    'AAG': 0, 'CTT': 0, 'AGA': 0, 'TCA': 0,
                    'TGA': 0, 'TCT': 0, 'CAT': 0, 'ATG': 0,
                    'AAC': 0, 'GTT': 0, 'ATC': 0, 'GAT': 0,
                    'ACA': 0, 'TGT': 0, 'AGT': 0, 'ACT': 0,
                    'CCA': 0, 'TGG': 0, 'GTA': 0, 'TAC': 0,
                    'TAG': 0, 'CTA': 0, 'GGA': 0, 'TCC': 0,
                    'CAG': 0, 'CTG': 0, 'GCA': 0, 'TGC': 0,
                    'ACC': 0, 'AGC': 0, 'GGT': 0, 'GCT': 0,
                    'AGG': 0, 'CCT': 0, 'GAG': 0, 'CTC': 0,
                    'CAC': 0, 'GTG': 0, 'GAC': 0, 'GTC': 0,
                    'CGA': 0, 'TCG': 0, 'ACG': 0, 'CGT': 0,
                    'GCC': 0, 'GGC': 0, 'CCC': 0, 'GGG': 0,
                    'CCG': 0, 'CGG': 0, 'CGC': 0, 'GCG': 0}

        nt4_list = []
        for i in ['A', 'T', 'C', 'G']:
            for j in ['A', 'T', 'C', 'G']:
                for m in ['A', 'T', 'C', 'G']:
                    for n in ['A', 'T', 'C', 'G']:
                        nt4_list.append(i + j + m + n)

        nt4_dict = dict(zip(nt4_list, [0] * 256))

        offset = 0
        for chromosome in self.data:
            for nt in chromosome[1]:
                nt_dict[nt] += 1
                offset += 1

            for index in range(len(chromosome[1]) - 2):
                nt2_dict[chromosome[1][index:index + 2]] += 1

        for dictionary in [nt_dict, nt2_dict]:
            sum_of_count = sum(dictionary.values())
            # print("the sum of {} is: {}".format(dictionary,sum(dictionary.values())),offset[[nt_dict,nt2_dict,nt3_dict].index(dictionary)])
            for nt, count in dictionary.items():
                dictionary[nt] = count / sum_of_count

        print("nt_ratio is:", nt_dict)
        print("the sum of nt is:", offset)
        print("nt2_ratio is: ", nt2_dict)
        print("the sum of 2nt is:", sum(nt2_dict.values()))

        return None

    def get_genome_dict_ratio(self, library, output_to_excel=None):

        '''
        get the ratio of each nucleotide or each two nucleotides in the genome
        '''

        length_of_lib = len(library)
        nt_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        nt2_dict = {'AA': 0, 'AT': 0, 'AC': 0, 'AG': 0, 'TA': 0, 'TT': 0, 'TC': 0, 'TG': 0, 'CA': 0, 'CT': 0, 'CC': 0,
                    'CG': 0, 'GA': 0, 'GT': 0, 'GC': 0, 'GG': 0}
        nt3_dict = {'AAA': 0, 'TTT': 0, 'AAT': 0, 'ATT': 0,
                    'ATA': 0, 'TAT': 0, 'GAA': 0, 'TTC': 0,
                    'CAA': 0, 'TTG': 0, 'TAA': 0, 'TTA': 0,
                    'AAG': 0, 'CTT': 0, 'AGA': 0, 'TCA': 0,
                    'TGA': 0, 'TCT': 0, 'CAT': 0, 'ATG': 0,
                    'AAC': 0, 'GTT': 0, 'ATC': 0, 'GAT': 0,
                    'ACA': 0, 'TGT': 0, 'AGT': 0, 'ACT': 0,
                    'CCA': 0, 'TGG': 0, 'GTA': 0, 'TAC': 0,
                    'TAG': 0, 'CTA': 0, 'GGA': 0, 'TCC': 0,
                    'CAG': 0, 'CTG': 0, 'GCA': 0, 'TGC': 0,
                    'ACC': 0, 'AGC': 0, 'GGT': 0, 'GCT': 0,
                    'AGG': 0, 'CCT': 0, 'GAG': 0, 'CTC': 0,
                    'CAC': 0, 'GTG': 0, 'GAC': 0, 'GTC': 0,
                    'CGA': 0, 'TCG': 0, 'ACG': 0, 'CGT': 0,
                    'GCC': 0, 'GGC': 0, 'CCC': 0, 'GGG': 0,
                    'CCG': 0, 'CGG': 0, 'CGC': 0, 'GCG': 0}

        nt4_list = []
        for i in ['A', 'T', 'C', 'G']:
            for j in ['A', 'T', 'C', 'G']:
                for m in ['A', 'T', 'C', 'G']:
                    for n in ['A', 'T', 'C', 'G']:
                        nt4_list.append(i + j + m + n)

        nt4_dict = dict(zip(nt4_list, [0] * 256))

        for base4 in library:
            nt4_dict[base4] += len(library[base4])

            for sequence in library[base4]:

                for nt in base4 + sequence:
                    nt_dict[nt] += 1

                for index_of_base3 in range(0, len(sequence), 3):
                    nt3_dict[sequence[index_of_base3:index_of_base3 + 3]] += 1

        print("the sum of library's keys is:", len(library.keys()))
        print("the sum of library's values is:", sum([len(i) for i in library.values()]))

        for index_of_dictionary, dictionary in enumerate([nt_dict, nt3_dict, nt4_dict]):
            sum_of_count = sum(dictionary.values())

            print("the count of {} is:{}".format(["nt_dict", "nt3_dict", "nt4_dict"][index_of_dictionary], dictionary))

            # print("the sum of {} is: {}".format(dictionary,sum(dictionary.values())),offset[[nt_dict,nt2_dict,nt3_dict].index(dictionary)])
            for nt, count in dictionary.items():
                dictionary[nt] = count / sum_of_count
            print("the ratio of {}  is:{}".format(["nt_dict", "nt3_dict", "nt4_dict"][index_of_dictionary], dictionary))
            # dictionary = dict(sorted(dictionary.items(),key=lambda d:d[1],reverse=1))

        print(nt3_dict, nt4_dict)
        nt3_dict = dict(sorted(nt3_dict.items(), key=lambda d: d[1], reverse=1))
        nt4_dict = dict(sorted(nt4_dict.items(), key=lambda d: d[1], reverse=1))
        print(nt_dict, nt3_dict, nt4_dict)
        return nt_dict, nt3_dict, nt4_dict

    def _loadData(self):
        '''
        load data from self.pth file

        No arguments and returns
        '''
        try:
            fin = open(self.pth, 'r')
        except:
            sys.stderr.write('Can\'t open genome file %s\n' % self.pth)
            exit(int(sys._getframe().f_lineno) + 20000)
        else:
            name = ''
            seq = ''
            for line in fin:
                line = line.strip('\n\r')
                m = re.search('^>(.+)$', line)
                if m:  # name
                    if seq != '':
                        self.data.append((name, seq))
                        name = ''
                        seq = ''
                    name = m.group(1)
                    continue
                line = line.upper()
                if not set(line).issubset(set('ATCG')):
                    line = re.sub(r'[^ATCG]', "", line)
                seq += line
            if seq != '':
                self.data.append((name, seq))
        finally:
            fin.close()
            if self.include_N:
                print("The genome include base N.")

        offset = 0
        for chromosome in self.data:
            offset += len(chromosome[1])

        self.length = offset
        print("The len of genome is:{}".format(offset))
        return

    def find(self, subseq):
        '''
        find index of subseq in genome

        Arguments:
            self --- object self pointer
            subseq --- subseq

        Returns:
            idx --- index of subseq in genome, if not find return -1
        '''
        if len(subseq) == self.libray_size:
            if subseq in self.library:
                idx = len(self.library[subseq]) - 1
                return self.library[subseq][random.randint(0, idx)]
            else:
                return -1
        else:
            idx = 0
            for chromosome in self.data:
                ix = chromosome[1].find(subseq)
                if ix == -1:
                    idx += len(chromosome[1])
                else:
                    return idx + ix
            return -1

    def __getitem__(self, item):
        '''
        get sub seq from genome by item
        item can be an index of a nucleotide or a slice

        Arguments:
            self --- object self pointer
            item --- index or a slice of index

        Return:
            subseq or nucleotide
        '''

        start = None
        end = None
        if isinstance(item, slice):
            start = item.start
            end = item.stop
            if end is None or start is None:
                sys.stderr.write('*** genome object refuse a slice with ')
                sys.stderr.write('no border ***\n')
                exit(20079)
            for chromosome in self.data:
                if start >= len(chromosome[1]):
                    start -= len(chromosome[1])
                    end -= len(chromosome[1])
                    continue
                else:
                    return chromosome[1][start:end]
        else:
            start = item
            for chromosome in self.data:
                if start >= len(chromosome[1]):
                    start -= len(chromosome[1])
                    continue
                else:
                    return chromosome[1][start]

        sys.stderr.write('*** Out of range when call ')
        sys.stderr.write('__getitem__(%s)***' % str(item))
        exit(20097)

    def construct_dict(self, size=None, seed=None, save_genome=False):
        genome_library_dict = {}
        gn = Genome(self.genome_path)
        for kmer in set(self.kmer_list):
            if size:
                lib = gn.construct_library(kmer, [size, seed])
                power = math.log(size, 4)
                print("power is {}".format(power), "library{} complete".format(kmer))
            else:
                lib = gn.construct_library(kmer)
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
                    print("\r{:.2%} has completed,this cycle cost: {:.2f}s".format(count / len(lib),
                                                                                   time.time() - cycle_init),
                          end='')
                    cycle_init = time.time()

            genome_nt_dict, genome_nt3_dict, genome_nt4_dict = gn.get_genome_dict_ratio(gn_dict)
            genome_data = [gn_dict, genome_nt_dict, genome_nt3_dict, genome_nt4_dict]
            print("gn.lib completed")
            # if not os.path.exists("G:\work_in_cngb\genome-size\genome_dict\{}".format(genome_name)):
            if save_genome:
                if not os.path.exists(self.library_output_path):
                    os.mkdir(self.library_output_path)
                if size:
                    np.save(os.path.join(self.library_output_path,
                                         "genome_lib-{}-seed_{}-power_{}-bioseq.npy".format(kmer, seed, power)),
                            genome_data)
                else:
                    np.save(os.path.join(self.library_output_path, "genome_lib-{}-bioseq.npy".format(kmer)),
                            genome_data)
            genome_library_dict[kmer] = genome_data
        return genome_library_dict