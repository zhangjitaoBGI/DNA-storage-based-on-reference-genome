# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:59:38 2021

@author: Jitao Zhang
"""
import numpy as np
import math
import struct
import random
import copy
import os

segment_length = 543
dest = r"E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/result/jpeg_compression_test/Div2k/length_{}/".format(segment_length)

# with open("../image/{}.jpg".format(image_name),'rb') as fn:
image_folder_path = r"E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\datasets\DIV2K_train_HR\small-128\jpeg_75"
image_list = os.listdir(image_folder_path)[:200]

for image in image_list:
    image_path = image_folder_path + "/" + image
    with open(image_path,'rb') as fn:
        fn1 = fn.read()
        size = len(fn1)
        row_count = math.ceil(size*8/segment_length)
        matrix = [[0 for _ in range(segment_length)] for _ in range(row_count)]
        supplementary_bit = segment_length - size*8%segment_length
        
        row = 0
        col = 0
        for byte_index in range(size):
            # Read a file as bytes
            # one_byte = file.read(1)
            one_byte = fn1[byte_index:byte_index + 1]
            element = list(map(int, list(str(bin(struct.unpack("B", one_byte)[0]))[2:].zfill(8))))
            for bit_index in range(8):
                matrix[row][col] = element[bit_index]
                col += 1
                if col == segment_length:
                    col = 0
                    row += 1
            
        for i in range(10):
            bit_sequence = ''
            byte_sequence= b''
            maxmium_row_of_file_header = int(int(0x026e)*8/segment_length)+1
            missing_row = np.random.choice(range(maxmium_row_of_file_header,row_count),math.ceil((row_count-1)*0.01),replace=False)
    
            print(i,missing_row)
            new_matrix = copy.deepcopy(matrix)
            for row in missing_row:
                random_sequence = [random.randint(0,1) for _ in range(segment_length)]
                new_matrix[row] = random_sequence
            for row in new_matrix:
                bit_sequence += ''.join([str(j) for j in row])
            
            bit_sequence = bit_sequence[:-supplementary_bit]
            for one_byte in range(0,len(bit_sequence),8):
                #print(bit_sequence[one_byte:one_byte+8],int(bit_sequence[one_byte:one_byte+8],2))
                byte_sequence += int(bit_sequence[one_byte:one_byte+8],2).to_bytes(1, byteorder='little', signed=False)
        
            dest_folder = dest+image[:-4]
            if not os.path.exists(dest_folder):
                os.makedirs(dest_folder)
            
            with open(dest_folder+"/{:0>4}.jpg".format(i),'wb') as fn_out:
                fn_out.write(byte_sequence)
        