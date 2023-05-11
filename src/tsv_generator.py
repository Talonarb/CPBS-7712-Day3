# This script takes the final sequence composition dictionary, dictionary of sequence 25mer alignments, sequence ID and Sequence dictionary, Position_dict, and Query as string
# and returns a TSV of which indexes the sequence aligns to the query, as well as the indexes of the overlapping portion of the sequence within the sequence. This function starts
# at the first matching set of sequence of 25, extends that until a mismatch occurs, and records the sequence_ID, name of longest contig, starting position of the matched sub-sequence,
# ending position of the subsequence, and the start and end of the same sub-sequence within the query. 
import pandas as pd
def comp_tsv_generator(Final_Seq_Comp_Dict, Seq_24mer_aligned, Seq_ID_Seq_Dict,Position_Dict, Fasta_String, suggested_bp_overlap_minimum):
    dict_for_tsv = {}
    contig_comp_list = []
    for i in range(len(Fasta_String)):        
        key1 = i
        key2_holder = list(Final_Seq_Comp_Dict[i])
        key2 = key2_holder[0]
        thing_to_get_seqs_from = list(Final_Seq_Comp_Dict[key1][key2])
        for i in thing_to_get_seqs_from:
            seqs = i
            if seqs[0] == '>':
                contig_comp_list.append(seqs)
    contig_comp_list = list(set(contig_comp_list))
    seq_ID_list = []
    qseq_id_list = []
    seq_start_list = []
    seq_end_list = []
    qstart_list = []
    qend_list = []
    for i in contig_comp_list:
        Sequence_ID = i
        Sequence_of_ID = str(Seq_ID_Seq_Dict[i])
        mer_aligned = str([Seq_24mer_aligned[Sequence_ID]])
        mer_aligned = mer_aligned[2:-2]
        position_of_mer_in_seq = Sequence_of_ID.find(mer_aligned)
        position_of_mer_in_query = Position_Dict[mer_aligned]
        counter = 0
        try:
            query_to_match = Fasta_String[position_of_mer_in_query + suggested_bp_overlap_minimum -1]
            seq_to_match = Sequence_of_ID[position_of_mer_in_seq + suggested_bp_overlap_minimum -1]
            while seq_to_match == query_to_match:
                try:
                    # The test portion is used to trigger an early index error and stop the counter from incrementing
                    seq_to_match_test = Sequence_of_ID[position_of_mer_in_seq + suggested_bp_overlap_minimum + counter + 1 ]
                    query_to_match_test = Fasta_String[position_of_mer_in_query + suggested_bp_overlap_minimum + counter + 1 ]
                    counter +=1
                    query_to_match = Fasta_String[position_of_mer_in_query + suggested_bp_overlap_minimum + counter]
                    seq_to_match = Sequence_of_ID[position_of_mer_in_seq + suggested_bp_overlap_minimum + counter]
                except IndexError:
                    break
            seq_ID_list.append(Sequence_ID)
            seq_start_list.append(position_of_mer_in_seq)
            seq_end_list.append(position_of_mer_in_seq + suggested_bp_overlap_minimum + counter)
            qseq_id_list.append('contig1')
            qstart_list.append(position_of_mer_in_query)
            qend_list.append(position_of_mer_in_query + counter + suggested_bp_overlap_minimum)
        except IndexError:
            seq_ID_list.append(Sequence_ID)
            seq_start_list.append(position_of_mer_in_seq)
            seq_end_list.append(position_of_mer_in_seq + suggested_bp_overlap_minimum + counter)
            qseq_id_list.append('contig1')
            qstart_list.append(position_of_mer_in_query)
            qend_list.append(position_of_mer_in_query + suggested_bp_overlap_minimum + counter)
            continue
    dict_for_tsv['sseqid'] = seq_ID_list
    dict_for_tsv['qseqid'] = qseq_id_list
    dict_for_tsv['sstart'] = seq_start_list
    dict_for_tsv['send'] = seq_end_list
    dict_for_tsv['qstart'] = qstart_list 
    dict_for_tsv['qend'] = qend_list
    dataframe_for_export = pd.DataFrame(dict_for_tsv)
    dataframe_for_export.to_csv('alleles.aln', sep="\t")

                


