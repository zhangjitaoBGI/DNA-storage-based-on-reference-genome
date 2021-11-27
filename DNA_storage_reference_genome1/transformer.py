import numpy as np
import genome
import copy
def get_nt_to_pixel_data_dict(encode_matrix_shape, nt3_ratio_list, diff_dict):
    diff_transforming_dict = {"Y_channel": {}, "Cb_channel": {}, "Cr_channel": {}}
    diff_dict_new =  copy.deepcopy(diff_dict)
    for channel in diff_dict_new:
        if len(diff_dict_new[channel]) != 64:
            print("length of diff_dict[{}] not equal to 64.".format(channel))
            for diff in range(31, -1, -1):

                if diff not in diff_dict_new[channel]:
                    diff_dict_new[channel].append(diff)

                if 63 - diff not in diff_dict_new[channel]:
                    diff_dict_new[channel].append(63 - diff)

    for index_of_channel, channel in enumerate(diff_transforming_dict):
        if encode_matrix_shape[index_of_channel] == 2:
            diff_transforming_dict[channel] = dict(zip(nt3_ratio_list[0].keys(), diff_dict_new[channel]))
        elif encode_matrix_shape[index_of_channel] == 4:
            diff_transforming_dict[channel] = dict(zip(nt3_ratio_list[1].keys(), diff_dict_new[channel]))
            
    print("diff_transforming_dict complete.")
    return diff_transforming_dict


def convert_genome_quadruplets(genome_dict):
    new_dict = {}
    for key in genome_dict:
        new_dict[genome.bioseq_convert_int(key)] = genome_dict[key]
    genome_dict = new_dict
    return genome_dict

def get_transformed_genome_dict(encode_matrix_shape, small_genome_dict, large_genome_dict, diff_transforming_dict,decode_mode=False,mean_pixel_dict=None):
    transformed_genome_dict = {"Y_channel": {}, "Cb_channel": {}, "Cr_channel": {}}
    if decode_mode:
        transformed_genome_dict_decode = {"Y_channel": [], "Cb_channel": [], "Cr_channel": []}
    small_genome_dict = convert_genome_quadruplets(small_genome_dict)
    large_genome_dict = convert_genome_quadruplets(large_genome_dict)

    for index_of_channel, channel in enumerate(transformed_genome_dict):
        if encode_matrix_shape[index_of_channel] == 2:
            genome_dict = small_genome_dict
        elif encode_matrix_shape[index_of_channel] == 4:
            genome_dict = large_genome_dict
        if mean_pixel_dict:
            for mean in mean_pixel_dict[channel]:
                if mean not in genome_dict:
                    genome_lib_keys = np.array(list(genome_dict.keys()))
                    candidate_mean_index_list = \
                        np.where(np.abs(genome_lib_keys - mean) == min(np.abs(genome_lib_keys - mean)))[0]
                    candidate_mean_list = [genome_lib_keys[index] for index in candidate_mean_index_list]
                    for mean in candidate_mean_list:
                        if mean not in transformed_genome_dict[channel]:
                            transformed_genome_dict[channel][mean] = []
                            for one_base3_unit in genome_dict[mean]:
                                one_base3_unit = [one_base3_unit[i:i + 3] for i in
                                                  range(0, (encode_matrix_shape[index_of_channel] ** 2) * 3, 3)]
                                base3_to_diff_unit = [diff_transforming_dict[channel][base3] for base3 in
                                                      one_base3_unit]
                                transformed_genome_dict[channel][mean].append(np.array(base3_to_diff_unit))
                else:
                    transformed_genome_dict[channel][mean] = []
                    for one_base3_unit in genome_dict[mean]:
                        one_base3_unit = [one_base3_unit[i:i + 3] for i in
                                          range(0, (encode_matrix_shape[index_of_channel] ** 2) * 3, 3)]
                        if "" in [base3 for base3 in one_base3_unit]:
                            print(channel,genome_dict[mean][0])
                            print(one_base3_unit)
                        base3_to_diff_unit = [diff_transforming_dict[channel][base3] for base3 in one_base3_unit]
                        transformed_genome_dict[channel][mean].append(np.array(base3_to_diff_unit))
        else:
            for mean in genome_dict:
                transformed_genome_dict[channel][mean] = []
                for one_base3_unit in genome_dict[mean]:
                    one_base3_unit = [one_base3_unit[i:i + 3] for i in
                                      range(0, (encode_matrix_shape[index_of_channel] ** 2) * 3, 3)]
                    base3_to_diff_unit = [diff_transforming_dict[channel][base3] for base3 in one_base3_unit]
                    transformed_genome_dict[channel][mean].append(np.array(base3_to_diff_unit))
            
    
        if decode_mode:
            for index_of_mean,mean in enumerate(transformed_genome_dict[channel]):
                if index_of_mean==0:
                    transformed_genome_dict_decode[channel] = np.int32(mean) + transformed_genome_dict[channel][mean]
                else:                    
                    transformed_genome_dict_decode[channel] = np.append(transformed_genome_dict_decode[channel],
                                                                        np.int32(mean) + transformed_genome_dict[channel][mean],axis=0)
            #transformed_genome_dict_decode[channel] = np.array(transformed_genome_dict_decode[channel])
    
    if decode_mode:
        print("transformed_genome_dict_decode complete.")        
        return transformed_genome_dict_decode
    else:
        print("transformed_genome_dict complete.")
        return transformed_genome_dict
