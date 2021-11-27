import sys
sys.path.append("../DNA_storage_YYC")
sys.path.append("../DNA_storage_YYC/utils")

import codec_factory
import yyc
import utils
def YYC_encoding(compressed_seq,dna_path,model_path,segment_length):
    [support_base, rule1, rule2] = ["A", [0, 1, 0, 1], [[1, 1, 0, 0], [1, 0, 0, 1], [1, 1, 0, 0], [1, 1, 0, 0]]]
    tool = yyc.YYC(support_bases=support_base, base_reference=rule1, current_code_matrix=rule2,
                   search_count=100, max_homopolymer=4, max_content=0.6)
    
    '''
    max_matrix_size = int(len(compressed_seq)*8/segment_length[0])
    
    if max_matrix_size%2 != 0:
        max_matrix_size -= 1
    compressed_seq_list = [compressed_seq[:int(max_matrix_size*segment_length[0]/8)],compressed_seq[int(max_matrix_size*segment_length[0]/8):]] 
    
    dna_sequences1 = codec_factory.encode(
        method=tool,
        index_seq_bytes=compressed_seq_list[0],
        output_path=dna_path[0],
        model_path=model_path[0],
        segment_length = segment_length[0],
        need_index=True,
        need_log=True
    )
    
    dna_sequences2 = codec_factory.encode(
    method=tool,
    index_seq_bytes=compressed_seq_list[1],
    output_path=dna_path[1],
    model_path=model_path[1],
    segment_length = segment_length[1],
    need_index=True,
    need_log=True
)    
    '''
    dna_sequences = codec_factory.encode(
    method=tool,
    index_seq_bytes=compressed_seq,
    output_path=dna_path,
    model_path=model_path,
    segment_length = segment_length,
    need_index=True,
    need_log=True
    )
if __name__ == '__main__':
    print("a")