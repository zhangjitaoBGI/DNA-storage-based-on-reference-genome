import sys
sys.path.append("..")
sys.path.append("../../DNA_storage_YYC")
import time
import genome_storage
import count
import encode_method1
import decode_method1
import transformer
import numpy as np

img_name = "lenna"
image_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/image/{}.png".format(img_name)
genome_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/fna/S288C.fna"
txid = 4932    #txid is the id number of genome in the NCBI genome database
txid_and_genome_dict = {txid:genome_path}
library_output_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/genome_dict/S288C"
encode_matrix_shape = [2,4,4]
dna_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/output/S288C/{}_random.dna".format(img_name)
model_path = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/output/S288C/{}_random.pkl".format(img_name)
segment_length = 543
dest = "E:/work in CNGB/DNA_storage/DNA_storage_reference_genome1/example1/result/{}_test_random.bmp".format(img_name)

if __name__ == '__main__':
    start_time = time.time()
    
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
    transformed_genome_dict = transformer.get_transformed_genome_dict(encode_matrix_shape, small_genome_dict, large_genome_dict, diff_transforming_dict)
    seq_in_genome, total_score_in_genome, seq_of_index = tool.encoded_seq_in_genome(rearranged_data_dict, transformed_genome_dict)
    compressed_seq = tool.connecting_file_header_to_index_seq(txid,diff_dict,seq_of_index,segment_length)
    np.save("compressed_seq.npy",compressed_seq)
    
    compressed_seq = np.load("compressed_seq_{}.npy".format(img_name),allow_pickle=True).item()
    encode_method1.YYC_encoding(compressed_seq,dna_path,model_path,segment_length=segment_length)
    del tool,rearranged_data_dict,diff_transforming_dict,transformed_genome_dict,seq_in_genome, total_score_in_genome, seq_of_index
    #decode_method.decode(dna_path, model_path, txid_and_genome_dict, library_output_path,dest,compressed_seq)
    decode_method1.decode(dna_path, model_path, txid_and_genome_dict, library_output_path,dest,compressed_seq)
    print("This image cost {} s".format(time.time() - start_time))
    