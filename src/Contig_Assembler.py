# Helper function to match list indices between a list of nucleotides and a list of sequence_IDs for tracking, based on https://datagy.io/python-list-find-all-index/
def index_tracker(values_as_list, most_frequent_nucleotide):
    indices = []
    for idx, value in enumerate(values_as_list):
        if value == most_frequent_nucleotide:
            indices.append(idx)
    return indices 

# Helper function to get a the largest list of consecutive numbers for contig assembly, via https://stackoverflow.com/questions/8664708/how-does-one-find-the-largest-consecutive-set-of-numbers-in-a-list-that-are-not

def max_subsequence_consecutive(position_list):
    from itertools import groupby
    from operator import itemgetter
    longest_list = []
    for k, g in groupby(enumerate(position_list), lambda ix : ix[0] - ix[1]):
        check_list = list(map(itemgetter(1), g))
        if len(check_list) > len(longest_list):
            longest_list = check_list
    return longest_list

# Function to assemble a contig with a dictionary of positions as keys and a list of both Sequences and the nucleotide they contain at a given position as values 
def contig_assembly(Contig_Dict, Relative_Position_of_Nucleotide_in_Sequence_dict):
    from statistics import mode
    position_list = []
    final_contig = ""
    final_seq_comp_dict = {}
    for i in Contig_Dict.keys():
        # Skips any positions with a coverage of less than 4, as the sequence ID and nucleotide are stored, divides by a factor of 8 instead of 4
        if (len((Contig_Dict[i]))) < 8:
            continue
        # Calculates the most represented nucleotide at a position, and returns the position and sequence_ID from which the position is represented
        else:
            values_as_list = list(Contig_Dict[i])
            Seq_ID_list = values_as_list[::2]
            nucleotide_list = values_as_list[1::2]
            most_frequent_nucleotide = mode(nucleotide_list)
            seq_ID_indexes = index_tracker(nucleotide_list, most_frequent_nucleotide)
        # Skips any positions with a coverage of less than 4 after checking the number of nucleotide occurences
            if len(seq_ID_indexes) < 4:
                continue
            else:
        # construct the list of sequence IDs corresponding to the most common nucleotide and create a list of all positions with coverage
                position_list.append(i)
                seq_IDs_for_contig = [Seq_ID_list[x] for x in seq_ID_indexes]
                for y in seq_IDs_for_contig:
                    final_seq_comp_dict[i] = {most_frequent_nucleotide:Relative_Position_of_Nucleotide_in_Sequence_dict[i][most_frequent_nucleotide]}
    position_list.sort(key=int)
    longest_contig_positions = max_subsequence_consecutive(position_list)
    for i in longest_contig_positions:
        values_as_list = list(Contig_Dict[i])
        Seq_ID_list = values_as_list[::2]
        nucleotide_list = values_as_list[1::2]
        most_frequent_nucleotide = mode(nucleotide_list)
        seq_ID_indexes = index_tracker(nucleotide_list, most_frequent_nucleotide)
        final_contig+= most_frequent_nucleotide
    # deleting the keys from the final_seq_comp_dict that are not in the longest contig position list, as there were issues making a dict after finding max subsequence
    for key in list(final_seq_comp_dict.keys()):
      if key not in longest_contig_positions:
         del final_seq_comp_dict[key]
    return final_seq_comp_dict, final_contig


        
    

            


