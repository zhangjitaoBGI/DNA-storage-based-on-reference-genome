def get_difference_distribution(rearranged_data_list):
    mean_pixel_dict = {'Y_channel': {}, 'Cb_channel': {}, 'Cr_channel': {}}
    diff_dict = {'Y_channel': {}, 'Cb_channel': {}, 'Cr_channel': {}}
    for index_of_channel, channel in enumerate(diff_dict):
        for index in range(len(rearranged_data_list[channel])):
            mean_pixel = rearranged_data_list[channel][index][0]
            if mean_pixel in mean_pixel_dict[channel]:
                mean_pixel_dict[channel][mean_pixel] += 1
            else:
                mean_pixel_dict[channel][mean_pixel] = 1

            for diff in rearranged_data_list[channel][index][1]:
                diff = int(diff)
                if diff in diff_dict[channel]:
                    diff_dict[channel][diff] += 1
                else:
                    diff_dict[channel][diff] = 1

        print("The sum of nt in {} channel is {}".format(channel, sum(diff_dict[channel].values())))


        for dictionary in [mean_pixel_dict, diff_dict]:
            sum_of_values = sum(dictionary[channel].values())

            for key, value in dictionary[channel].items():
                dictionary[channel][key] = value / sum_of_values

            dictionary[channel] = dict(sorted(dictionary[channel].items(), key=lambda d:d[1], reverse=1))
            dictionary[channel] = list(dictionary[channel].keys())
        '''
        sum_of_values = sum(diff_dict[channel].values())

        for diff,count in diff_dict[channel].items():
            diff_dict[channel][diff] = count/sum_of_values

        diff_dict[channel] = dict(sorted(diff_dict[channel].items(),key=lambda d:d[0]))
        '''
    return mean_pixel_dict, diff_dict