import sys
sys.path.append("..")
sys.path.append(".")
sys.path.append("../DNA_storage_YYC")
import codec_factory
import genome_storage
import numpy as np
import cv2
import transformer
import random
import copy

def YYC_decoding(dna_path,model_path,loss_seq_ratio=None):
    dna_path = dna_path
    model_path = model_path
    _, decoded_bytes = codec_factory.decode(
        model_path=model_path,
        input_path=dna_path,
        has_index=True,
        need_log=True,
        need_out_put_index=True,
        loss_seq_ratio=loss_seq_ratio
    )

    return decoded_bytes


def get_header_and_index(decoded_bytes):
    decoded_index = []  
    header_length = int.from_bytes(decoded_bytes[:2],byteorder="little", signed=False)
    #segment_length[0] = int.from_bytes(decoded_bytes[2:4],byteorder="little", signed=False)
    #segment_length[1] = int.from_bytes(decoded_bytes[4:6],byteorder="little", signed=False)
    txid = int.from_bytes(decoded_bytes[2:6], byteorder="little", signed=False)
    image_size = [int.from_bytes(decoded_bytes[i:i+2], byteorder="little", signed=False) for i in [6,8]]
    encode_matrix = [int.from_bytes(decoded_bytes[i:i+1], byteorder="little", signed=False) for i in [10,11,12]]
    diff_ratio_count = [int.from_bytes(decoded_bytes[i:i+1], byteorder="little", signed=False) for i in [13,14,15]]
    diff_dict = {"Y_channel":[],"Cb_channel":[],"Cr_channel":[]}
    off_set = 16
    for index_of_channel,channel in enumerate(diff_dict):
        diff_dict[channel] = [int.from_bytes(decoded_bytes[i:i+1], byteorder="little", signed=False) for i in range(off_set,off_set+diff_ratio_count[index_of_channel])]
        off_set += diff_ratio_count[index_of_channel]
    if off_set != header_length:
        print("header is error !")

    index_seq = decoded_bytes[header_length:]

    for index in range(0, len(index_seq), 3):
        decoded_index.append(
            int.from_bytes(index_seq[index:index + 3], byteorder="little",signed=False))
    return txid,image_size,encode_matrix,diff_dict,decoded_index


def retrieve_index_in_genome_library(decoded_index,image_size,encode_matrix,transformed_genome_dict):
    decoded_index_dict = {'Y_channel': [], 'Cb_channel': [], 'Cr_channel': []}
    one_channel_data_length = [int(image_size[0] * image_size[1] / (i ** 2)) for i in
                               encode_matrix]
    index_count_dict = {'Y_channel': decoded_index[:one_channel_data_length[0]],
                           'Cb_channel': decoded_index[one_channel_data_length[0]:one_channel_data_length[0] +
                                                                                   one_channel_data_length[1]],
                           'Cr_channel': decoded_index[one_channel_data_length[0] + one_channel_data_length[1]:]}
    for channel in index_count_dict:
        new_index_list = copy.deepcopy(index_count_dict[channel])
        for index_of_index,index in enumerate(index_count_dict[channel]):
            if index >= len(transformed_genome_dict[channel]):
                if index_of_index != 0:                    
                    for j in range(1,index_of_index+1):
                        if new_index_list[index_of_index-j] < len(transformed_genome_dict[channel]):
                            index = new_index_list[index_of_index-j]
                            new_index_list[index_of_index] = index
                            break
                        # elif j==(index_of_index-1):
                        #     print("There is no correct index in genome index.")
                        #     print(index_count_dict[channel])
                else:
                    for j in range(1,len(transformed_genome_dict[channel])):
                        if index_count_dict[channel][index_of_index+j] < len(transformed_genome_dict[channel]):
                            index = index_count_dict[channel][index_of_index+j]
                            new_index_list[index_of_index] = index
                            break
                    #index = random.randint(0,len(transformed_genome_dict[channel]))    #纯随机替代序列
            if index >= len(transformed_genome_dict[channel]):
                print("index is : ",index)
                print("index_of_index is : ",index_of_index)
                print("index_count_dict[channel] is : ", index_count_dict[channel][:index_of_index],"new_index_list is : ",new_index_list)
                print("index is not in transformed_genome_dict[channel]")
            decoded_index_dict[channel].append(transformed_genome_dict[channel][index])
            
    return decoded_index_dict


def decode_into_image(decoded_index_dict, dest,image_size,encode_matrix):
    decoded_file_data_dict = {'Y_channel': [], 'Cb_channel': [], 'Cr_channel': []}

    # length_of_block_in_3_channel = sum([int(i**2*16/i**2) for i in self.encode_matrix_shape])

    # for index in range(0,len(compressed_seq),length_of_block_in_3_channel):

    for index_of_channel, channel in enumerate(decoded_file_data_dict):
        #matrix_size = encode_matrix[index_of_channel] ** 2
        for index_of_pixel in range(0, len(decoded_index_dict[channel])):
            #mean_pixel = np.int16(decoded_index_dict[channel][index_of_pixel])
            if type(decoded_index_dict[channel][index_of_pixel]) is not np.ndarray:
                decoded_index_dict[channel][index_of_pixel] = np.array(decoded_index_dict[channel][index_of_pixel],dtype="uint16")
            image_data_list = list(decoded_index_dict[channel][index_of_pixel] - 31)


            for data in image_data_list:
                if data < 0:
                    image_data_list[image_data_list.index(data)] = 0
                elif data > 255:
                    image_data_list[image_data_list.index(data)] = 255

            decoded_file_data_dict[channel] += image_data_list
    # print(decoded_file_data_list)
    data_of_image_format_list = {'Y_channel': [], 'Cb_channel': [], 'Cr_channel': []}

    for index_in_channel, channel in enumerate(decoded_file_data_dict):
        shape = encode_matrix[index_in_channel]
        for block in range(0, len(decoded_file_data_dict[channel]), image_size[0] * shape):
            for row in range(shape):
                for step in range(0, image_size[0] * shape, shape ** 2):
                    for column in range(shape):
                        data_of_image_format_list[channel].append(
                            decoded_file_data_dict[channel][block + row * shape + step + column])

        data_of_image_format_list[channel] = np.array(data_of_image_format_list[channel], dtype="uint8").reshape(
            image_size)

    img_ycc = cv2.merge([data_of_image_format_list['Y_channel'], data_of_image_format_list['Cr_channel'],
                         data_of_image_format_list['Cb_channel']])
    # img_ycc = img_ycc.astype('uint8')
    # a = img_ycc1 - img_ycc
    img = cv2.cvtColor(img_ycc, cv2.COLOR_YCR_CB2BGR)
    cv2.imwrite(dest, img)

    

def decode(dna_path,model_path,txid_and_genome_dict,library_output_path,dest,compressed_seq=None,loss_seq_ratio=None,simulation_count=None):
    '''
    decoded_bytes = b''
    decoded_bytes += YYC_decoding(dna_path[0],model_path[0])
    decoded_bytes += YYC_decoding(dna_path[1],model_path[1])
    '''
    if simulation_count != None:
        decoded_index_list = []
        for i in range(simulation_count,1000):
            decoded_bytes = YYC_decoding(dna_path,model_path,loss_seq_ratio)
            if compressed_seq:
                print("source digital file == target digital file: " + str(compressed_seq == decoded_bytes))            
            decoded_index_list.append(get_header_and_index(decoded_bytes)[-1])
    txid, image_size, encode_matrix_shape, diff_dict,_ = get_header_and_index(decoded_bytes)    
    genome_path = txid_and_genome_dict[txid]
    tool = genome_storage.Genome_storage(encode_matrix_shape,genome_path,library_output_path)
    genome_library_dict = tool.construct_dict(save_genome=True)
    nt_ratio_list = [0, 0]
    nt3_ratio_list = [0, 0]
    nt4_ratio_list = [0, 0]
    small_genome_dict,nt_ratio_list[0],nt3_ratio_list[0],nt4_ratio_list[0] = genome_library_dict[16]
    large_genome_dict,nt_ratio_list[1],nt3_ratio_list[1],nt4_ratio_list[1] = genome_library_dict[52]
    del genome_library_dict
    diff_transforming_dict = transformer.get_nt_to_pixel_data_dict(encode_matrix_shape, nt3_ratio_list, diff_dict)
    transformed_genome_dict = transformer.get_transformed_genome_dict(encode_matrix_shape, small_genome_dict, large_genome_dict, diff_transforming_dict,decode_mode=True,mean_pixel_dict=None)
    if decoded_index_list:
        for index,decoded_index in enumerate(decoded_index_list):
            decoded_index_dict = retrieve_index_in_genome_library(decoded_index,image_size,encode_matrix_shape,transformed_genome_dict)
            path = dest + "{:0>4}.bmp".format(index+simulation_count)
            decode_into_image(decoded_index_dict,path,image_size,encode_matrix_shape)
            print("image {:0>4} completed".format(index+simulation_count),end=" ")
    else:
        decoded_index_dict = retrieve_index_in_genome_library(decoded_index,image_size,encode_matrix_shape,transformed_genome_dict)
        decode_into_image(decoded_index_dict,path,image_size,encode_matrix_shape)
    
if __name__ == '__main__':
    img_name = "flowers"
    image_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/image/{}.png".format(img_name)
    genome_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/fna/S288C.fna"
    txid = 4932    #txid is the id number of genome in the NCBI genome database
    txid_and_genome_dict = {txid:genome_path}
    library_output_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/genome_dict/S288C"
    encode_matrix_shape = [2,4,4]
    dna_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/output/S288C/{}_original_yyc.dna".format(img_name)
    #dna_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/validation/{}.dna".format(img_name)
    model_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/output/S288C/{}_original_yyc.pkl".format(img_name)
    segment_length = 543
    dest = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/result/{}_original_yyc.bmp".format(img_name)
    #decode(dna_path, model_path, txid_and_genome_dict, library_output_path,dest)
    decoded_bytes = YYC_decoding(dna_path,model_path)
    # if compressed_seq:
    #     print("source digital file == target digital file: " + str(compressed_seq == decoded_bytes))
    txid, image_size, encode_matrix_shape, diff_dict, decoded_index = get_header_and_index(decoded_bytes)
    '''
    genome_path = txid_and_genome_dict[txid]
    tool = genome_storage.Genome_storage(encode_matrix_shape,genome_path,library_output_path)
    genome_library_dict = tool.construct_dict(save_genome=True)
    nt_ratio_list = [0, 0]
    nt3_ratio_list = [0, 0]
    nt4_ratio_list = [0, 0]
    small_genome_dict,nt_ratio_list[0],nt3_ratio_list[0],nt4_ratio_list[0] = genome_library_dict[16]
    large_genome_dict,nt_ratio_list[1],nt3_ratio_list[1],nt4_ratio_list[1] = genome_library_dict[52]
    del genome_library_dict
    diff_transforming_dict = transformer.get_nt_to_pixel_data_dict(encode_matrix_shape, nt3_ratio_list, diff_dict)
    transformed_genome_dict = transformer.get_transformed_genome_dict(encode_matrix_shape, small_genome_dict, large_genome_dict, diff_transforming_dict,decode_mode=True,mean_pixel_dict=None)
    '''
    transformed_genome_dict = np.load("transformed_genome_dict.npy",allow_pickle=True).item()
    decoded_index_dict = retrieve_index_in_genome_library(decoded_index,image_size,encode_matrix_shape,transformed_genome_dict)
    decode_into_image(decoded_index_dict,dest,image_size,encode_matrix_shape)