# -*- coding: utf-8 -*-
"""
Created on Sun Jul 25 19:04:53 2021

@author: Jitao Zhang
"""
file_name = "flowers"

original_file = r"E:\work in CNGB\DNA_storage\DNA_storage_reference_genome1\example1\output\S288C\{}_random.dna".format(file_name)
sequenced_file = r"E:\work in CNGB\DNA_storage\DNA_storage_reference_genome1\example1\output\S288C\{}_test.dna".format(file_name)

fn1 = open(original_file)
fn2 = open(sequenced_file)

img1 = fn1.read()
img2 = fn2.read()

if img1 == img2:
    print('1')
else:
    print('2')