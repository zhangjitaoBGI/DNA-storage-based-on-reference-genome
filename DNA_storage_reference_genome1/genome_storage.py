import sys
sys.path.append("./")
import cv2
import numpy as np
import genome
import math
import time
import os


class Genome_storage(object):
    def __init__(self, encode_matrix_shape,genome_path,library_output_path=None,image_path=None):
        """
        Initialize Genome_storage
        :param inform_count: count of inform.
        :param hash_function_type: type of hash function: BKDR, DJB, JS, RS
        """
        if image_path != None:
            self.data_matrix = {'Y_channel': [], 'Cb_channel': [], 'Cr_channel': []}
            image = cv2.imread(image_path)
            image_ycc = cv2.cvtColor(image, cv2.COLOR_BGR2YCR_CB)
            self.data_matrix['Y_channel'] = np.ndarray.flatten(image_ycc[:, :, 0].astype("int16"))
            self.data_matrix['Cb_channel'] = np.ndarray.flatten(image_ycc[:, :, 2].astype("int16"))
            self.data_matrix['Cr_channel'] = np.ndarray.flatten(image_ycc[:, :, 1].astype("int16"))
            print(self.data_matrix['Y_channel'].itemsize)
            self.image_size = image.shape[0:2]
        self.encode_matrix_shape = encode_matrix_shape
        self.kmer_list = [4+i**2*3 for i in encode_matrix_shape]
        self.genome_path = genome_path
        self.library_output_path = library_output_path


    def image_data_convert_to_seq(self):
        reshaped_data_matrix = {'Y_channel': [], 'Cb_channel': [], 'Cr_channel': []}
        row = 0
        for index_of_channel, channel in enumerate(self.data_matrix):
            index = 0
            while (index < len(self.data_matrix[channel])):
                if (index < (row + self.encode_matrix_shape[index_of_channel]) * self.image_size[1]):
                    for pixel in range(index, index + self.image_size[1], self.encode_matrix_shape[index_of_channel]):

                        data_segment_matrix = []
                        for row_in_matrix in range(self.encode_matrix_shape[index_of_channel]):
                            for column_in_matrix in range(self.encode_matrix_shape[index_of_channel]):
                                data_segment_matrix.append(self.data_matrix[channel][
                                                               pixel + row_in_matrix * self.image_size[
                                                                   1] + column_in_matrix])

                        mean_in_data_segment_matrix = np.int16(np.mean(data_segment_matrix))

                        reshaped_data_matrix[channel].append(mean_in_data_segment_matrix)
                        for element in data_segment_matrix:
                            if element - mean_in_data_segment_matrix < -31:
                                reshaped_data_matrix[channel].append(-31)
                            elif element - mean_in_data_segment_matrix > 32:
                                reshaped_data_matrix[channel].append(32)
                            else:
                                reshaped_data_matrix[channel].append(element - mean_in_data_segment_matrix)

                    index += self.encode_matrix_shape[index_of_channel] * self.image_size[1]

                else:
                    row += self.encode_matrix_shape[index_of_channel]

        print("reshaped_data_matrix length = ", len([[] + i for i in reshaped_data_matrix.values()]))


        rearranged_data_dict = {'Y_channel': [], 'Cb_channel': [], 'Cr_channel': []}

        for index_of_channel, channel in enumerate(reshaped_data_matrix):

            for index_data_matrix in range(0, len(reshaped_data_matrix[channel]), self.encode_matrix_shape[index_of_channel] ** 2 + 1):  # One mean and four difference
                # block = [LZW.compress_short.int_convert_bioseq(reshaped_data_matrix[channel][index_data_matrix],8)]
                block = [reshaped_data_matrix[channel][index_data_matrix]]
                # block.append(np.array(data_matrix[channel][index_data_matrix+1:index_data_matrix+self.encode_matrix_shape[index_of_channel]**2+1])+63)

                block.append(np.array(reshaped_data_matrix[channel][
                                      index_data_matrix + 1:index_data_matrix + self.encode_matrix_shape[
                                          index_of_channel] ** 2 + 1]) + 31)
                rearranged_data_dict[channel].append(block)

            print("blocks of " + channel + ":", rearranged_data_dict[channel][:10],
                  "\ncounts of blocks in " + channel + ":", len(rearranged_data_dict[channel]))

        return rearranged_data_dict
    
    
    def construct_dict(self,size=None,seed=None,save_genome=False):

        genome_library_dict = {}
        gn = genome.Genome(self.genome_path)
        for kmer in set(self.kmer_list):
            if size:
                power = math.log(size,4)
                if not os.path.exists(os.path.join(self.library_output_path,"genome_lib-{}-seed_{}-power_{}-bioseq.npy".format(kmer, seed, power))):
                    lib = gn.construct_library(kmer, [size, seed])
                    print("power is {}".format(power), "library{} complete".format(kmer))
                else:
                    genome_data = np.load(os.path.join(self.library_output_path,"genome_lib-{}-seed_{}-power_{}-bioseq.npy".format(kmer, seed, power)),allow_pickle=True)
                    genome_library_dict[kmer] = genome_data
                    print("This genome exists and it has been read.")
                    continue
            else:
                if not os.path.exists(os.path.join(self.library_output_path,"genome_lib-{}-bioseq.npy".format(kmer))):
                    lib = gn.construct_library(kmer)
                    print("the lenght of gn.library is :" + " " + str(len(lib)))
                else:
                    genome_data = np.load(os.path.join(self.library_output_path,"genome_lib-{}-bioseq.npy".format(kmer)),allow_pickle=True)
                    genome_library_dict[kmer] = genome_data
                    print("This genome exists and it has been read.")
                    continue
            
            
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

            genome_nt_dict, genome_nt3_dict, genome_nt4_dict = gn.get_genome_dict_ratio(gn_dict)
            genome_data = [gn_dict, genome_nt_dict, genome_nt3_dict, genome_nt4_dict]
            print("gn.lib completed")
            #if not os.path.exists("G:\work_in_cngb\genome-size\genome_dict\{}".format(genome_name)):
            if save_genome:
                if not os.path.exists(self.library_output_path):
                    os.mkdir(self.library_output_path)
                if size:
                    np.save(os.path.join(self.library_output_path,"genome_lib-{}-seed_{}-power_{}-bioseq.npy".format(kmer, seed, power)), genome_data)
                else:
                    np.save(os.path.join(self.library_output_path,"genome_lib-{}-bioseq.npy".format(kmer)), genome_data)
            genome_library_dict[kmer] = genome_data
        return genome_library_dict
        


    def encoded_seq_in_genome(self, rearranged_data_list, transformed_genome_dict):
        seq_in_genome = {'Y_channel': [], 'Cb_channel': [], 'Cr_channel': []}
        total_score_in_genome = {'Y_channel': 0, 'Cb_channel': 0, 'Cr_channel': 0}
        seq_of_index = {'Y_channel': [], 'Cb_channel': [], 'Cr_channel': []}

        for index_of_channel, channel in enumerate(rearranged_data_list):
            value_size_in_gn_library = {}
            offset = 0
            for key in transformed_genome_dict[channel]:
                value_size_in_gn_library[key] = offset
                offset += len(transformed_genome_dict[channel][key])

            for index_of_block in range(0, len(rearranged_data_list[channel])):
                max_score = 0
                dif_array = np.array([])
                min_dif = 0
                index_in_gn_lib = []
                mean = rearranged_data_list[channel][index_of_block][0]

                if mean not in transformed_genome_dict[channel]:
                    genome_lib_keys = np.array(list(transformed_genome_dict[channel].keys()))

                    candidate_mean_index_list = \
                    np.where(np.abs(genome_lib_keys - mean) == min(np.abs(genome_lib_keys - mean)))[0]
                    candidate_mean_list = [genome_lib_keys[index] for index in candidate_mean_index_list]
                    min_dif_list = [0, 0]
                    dif_array_list = [0, 0]
                    for index_of_candidate_mean, candidate_mean in enumerate(candidate_mean_list):
                        target_diff_seq = rearranged_data_list[channel][index_of_block][1] - (candidate_mean - mean)

                        dif_array_list[index_of_candidate_mean] = np.sum(
                            np.abs(transformed_genome_dict[channel][candidate_mean] - target_diff_seq), axis=1)
                        min_dif_list[index_of_candidate_mean] = np.min(dif_array_list[index_of_candidate_mean])

                    if min_dif_list[0] <= min_dif_list[1]:
                        mean = candidate_mean_list[0]
                        min_dif = min_dif_list[0]
                        dif_array = dif_array_list[0]

                    else:
                        mean = candidate_mean_list[1]
                        min_dif = min_dif_list[1]
                        dif_array = dif_array_list[1]

                else:
                    dif_array = np.sum(np.abs(transformed_genome_dict[channel][mean]
                               - rearranged_data_list[channel][index_of_block][1]), axis=1)

                    min_dif = np.min(dif_array)

                max_score = 63 * (self.encode_matrix_shape[index_of_channel] ** 2) - min_dif

                total_score_in_genome[channel] += max_score

                index_in_gn_lib = np.where(dif_array == min_dif)[0][0]

                seq_of_index[channel].append(value_size_in_gn_library[mean] + index_in_gn_lib)

                candidate_list = [mean] + list(transformed_genome_dict[channel][mean][index_in_gn_lib])

                seq_in_genome[channel].append(candidate_list)

                print("\rchannel is: {0}, index_of_block is: {1},max_score is: {2}, seq_with_max_score is: {3}.".format(
                    channel, index_of_block, max_score, candidate_list), end="")

            print("\rOne channel complete.",end='')

            total_score_in_genome[channel] /= 63 * self.image_size[0] * self.image_size[1]
        print("total_score_in_genome is:", total_score_in_genome)

        return seq_in_genome, total_score_in_genome, seq_of_index


    def connecting_file_header_to_index_seq(self, txid, diff_ratio,seq_of_index,segment_length,out_put_seq_length=False):
        compressed_seq = []
        data_header = b""
        index_seq = b""
        data_header += txid.to_bytes(4,byteorder='little', signed=False)
        data_header += self.image_size[0].to_bytes(2,byteorder='little', signed=False)
        data_header += self.image_size[1].to_bytes(2,byteorder='little', signed=False)

        for i in self.encode_matrix_shape:
            data_header += i.to_bytes(1, byteorder='little', signed=False)

        diff_ratio_count = [len(diff_ratio[channel]) for channel in diff_ratio]
        for ratio_count in diff_ratio_count:
            data_header += ratio_count.to_bytes(1, byteorder='little', signed=False)

        for channel in diff_ratio:
            for diff in diff_ratio[channel]:
                data_header += diff.to_bytes(1,byteorder='little', signed=False)

        header_length = len(data_header) + 2
        print("header_length is:",header_length)
        print("length of data_header is:", len(data_header))
        data_header = header_length.to_bytes(2,byteorder='little', signed=False) + data_header
        print("length of data_header is:", len(data_header))
        # number_of_block = [int(len(seq_in_genome[channel])/(self.encode_matrix_shape[index_of_channel])**2*3+4) for index_of_channel,channel in enumerate(seq_in_genome)]
        

        for index_of_channel, channel in enumerate(seq_of_index):
            size_of_each_index = len(seq_of_index[channel][0].tobytes())
            if size_of_each_index != 3:
                for index in seq_of_index[channel]:
                    index_seq += index.tobytes()[:3]
            else:
                index_seq += np.array(seq_of_index[channel]).flatten().tobytes()


        print("length of index_seq_bytes is:", len(index_seq))

        compressed_seq = data_header + index_seq
        # with open(write_file_path, "wb") as file:
        #     print("Write storage data to file: " + write_file_path)
        #
        #     for byte in range(len(compressed_seq)):
        #         file.write(byte)
        
        if out_put_seq_length:
            
            return compressed_seq,header_length

        else:            

            return compressed_seq




if __name__ == '__main__':
    print("a")