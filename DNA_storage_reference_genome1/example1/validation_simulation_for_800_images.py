import sys
sys.path.append(".")
sys.path.append("../")
sys.path.append("../../DNA_storage_YYC")
import time
import genome_storage
import count
import encode_method1
import decode_method1
import transformer
import numpy as np
import math
import copy
import random
from decode_method1 import retrieve_index_in_genome_library,decode_into_image
import os
import cv2
import struct

# img_name = "baboon"
# image_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/image/{}.png".format(img_name)
genome_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/fna/S288C.fna"
txid = 4932    #txid is the id number of genome in the NCBI genome database
txid_and_genome_dict = {txid:genome_path}
library_output_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/genome_dict/S288C"
encode_matrix_shape = [2,4,4]
# dna_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/output/S288C/{}_random.dna".format(img_name)
# dna_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/validation/{}.dna".format(img_name)
# model_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/output/S288C/{}_random.pkl".format(img_name)
segment_length = 543
# segment_length = 200

def generate_random_index(compressed_index_seq,each_index_length,segment_length,loss_seq_ratio,test_count):
    size = len(compressed_index_seq)

    # Set init storage matrix
    matrix = [[0 for _ in range(segment_length)] for _ in range(math.ceil(size * 8 / segment_length))]
    # matrix = np.zeros((math.ceil(size * 8 / segment_length),segment_length),dtype="uint16")
    size = len(compressed_index_seq)
    supplementary_bit = segment_length - size*8%segment_length
    index_sequence_list = []
    row = 0
    col = 0
    
    for byte_index in range(size):
        # Read a file as bytes
        # one_byte = file.read(1)
        one_byte = compressed_index_seq[byte_index:byte_index + 1]
        element = list(map(int, list(str(bin(struct.unpack("B", one_byte)[0]))[2:].zfill(8))))
        for bit_index in range(8):
            matrix[row][col] = element[bit_index]
            col += 1
            if col == segment_length:
                col = 0
                row += 1
    
    for i in range(test_count):
        bit_sequence = ''
        byte_sequence = b''
        index_list = []
        missing_row = np.random.choice(range(len(matrix)),math.ceil(len(matrix)*0.01),replace=False)
        
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
            byte_sequence += int(bit_sequence[one_byte:one_byte+8],2).to_bytes(1, byteorder="little")
        
        for index_of_one_index in range(0,len(byte_sequence),each_index_length):
            index_list.append(int.from_bytes(byte_sequence[index_of_one_index:index_of_one_index+each_index_length], byteorder="little"))
        index_sequence_list.append(index_list)
        
    return index_sequence_list

def image_retrieval(index_sequence_list,image_size,encode_matrix_shape,transformed_genome_dict,dest):
    transformed_genome_dict_decode = {"Y_channel":0,"Cb_channel":0,"Cr_channel":0}
    for channel in transformed_genome_dict:
        for index_of_mean,mean in enumerate(transformed_genome_dict[channel]):
            if index_of_mean==0:
                transformed_genome_dict_decode[channel] = np.int32(mean) + transformed_genome_dict[channel][mean]
            else:                    
                transformed_genome_dict_decode[channel] = np.append(transformed_genome_dict_decode[channel],
                                                                    np.int32(mean) + transformed_genome_dict[channel][mean],axis=0)    

    for index,decoded_index in enumerate(index_sequence_list):
        decoded_index_dict = retrieve_index_in_genome_library(decoded_index,image_size,encode_matrix_shape,transformed_genome_dict_decode)
        path = dest + "{:0>4}.png".format(index)
        decode_into_image(decoded_index_dict,path,image_size,encode_matrix_shape)
    
        print("image {:0>4} completed".format(index),end=" ")


if __name__ == '__main__':
    start_time = time.time()
    image_filefolder_path = r"E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\datasets\DIV2K_train_HR\small-128\color"
    image_list = os.listdir(image_filefolder_path)[150:200]
    dest_path = r"E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/result/validation_replace_simulation/Div2k_test/length_543"
    for image in image_list:
        image_path = image_filefolder_path + "\\" + image
        image_size = cv2.imread(image_path).shape[:2]
        tool = genome_storage.Genome_storage(encode_matrix_shape,genome_path,library_output_path,image_path)
        genome_library_dict = tool.construct_dict(save_genome=True)
        kmer_list = [4+i**2*3 for i in set(encode_matrix_shape)]
        nt_ratio_list = [0, 0]
        nt3_ratio_list = [0, 0] 
        nt4_ratio_list = [0, 0]
        small_genome_dict,nt_ratio_list[0],nt3_ratio_list[0],nt4_ratio_list[0] = genome_library_dict[min(kmer_list)]
        large_genome_dict,nt_ratio_list[1],nt3_ratio_list[1],nt4_ratio_list[1] = genome_library_dict[max(kmer_list)]
        del genome_library_dict
        
        rearranged_data_dict = tool.image_data_convert_to_seq()
        mean_pixel_dict, diff_dict = count.get_difference_distribution(rearranged_data_dict)
        diff_transforming_dict = transformer.get_nt_to_pixel_data_dict(encode_matrix_shape, nt3_ratio_list, diff_dict)
        time_to_start_transform = time.time()
        transformed_genome_dict = transformer.get_transformed_genome_dict(encode_matrix_shape, small_genome_dict, large_genome_dict, diff_transforming_dict)
        time_cost_to_transform = time.time() - time_to_start_transform
        print("transformed_genome_dict of {} has completed and it cost {} min".format(image,time_cost_to_transform/60))
        
        seq_in_genome, total_score_in_genome, seq_of_index = tool.encoded_seq_in_genome(rearranged_data_dict, transformed_genome_dict)
        compressed_seq,header_length = tool.connecting_file_header_to_index_seq(txid,diff_dict,seq_of_index,segment_length,out_put_seq_length=True)
        
        compressed_index_seq = compressed_seq[header_length:]
        index_sequence_list =  generate_random_index(compressed_index_seq,3,segment_length,0.01,10)     
        dest = dest_path + "/{}/".format(image[:-4])
        if not os.path.exists(dest):
            os.mkdir(dest)
        
        image_retrieval(index_sequence_list,image_size,encode_matrix_shape,transformed_genome_dict,dest)
        print("{} has completed, and it cost {} min.".format(image,(time.time() - start_time)/60))
        del transformed_genome_dict,compressed_seq,compressed_index_seq,index_sequence_list
        start_time = time.time() 
        
    
    '''                                                                        
    np.save("compressed_seq.npy",compressed_seq)
    compressed_seq = np.load("compressed_seq_{}.npy".format(img_name),allow_pickle=True).item()
    encode_method1.YYC_encoding(compressed_seq,dna_path,model_path,segment_length=segment_length)
    #del tool,rearranged_data_dict,diff_transforming_dict,transformed_genome_dict,seq_in_genome, total_score_in_genome, seq_of_index
    #decode_method1.decode(dna_path, model_path, txid_and_genome_dict, library_output_path,dest,compressed_seq)
    dest = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/result/validation_replace_simulation/{}/".format(img_name)
    decode_method1.decode(dna_path, model_path, txid_and_genome_dict, library_output_path,dest,loss_seq_ratio=0.01,simulation_count=0)
    print("This image cost {} s".format(time.time() - start_time))
    '''