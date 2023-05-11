# Function to create a dictionary of all positions relative to Query as base, with negative value keys indicating positions preceding the query and keys greater than query length
# indicating positions after the end of the query, with values indicating the source of a nucleotide as Seq_ID and the nucleotide. Counter acts as a position tracker for each
# nucleotide relative to the position it holds relative to the query sequence. Another counter tracks the position of a nucleotide within a sequence and stores the position as key,
# with Seq_ID and position of a nucleotide within a sequence as value.
def Contig_Dict_Maker(forward_seqs_alignment_dict, Seq_ID_Seq_Dict, Position_Dict):  
	Contig_Dict = {}
	relative_position_of_nucleotide_in_seq_dict = {}
	for key in Seq_ID_Seq_Dict:
		position_of_match = forward_seqs_alignment_dict[key]
		inverted_position_dict = {v: k for k, v in Position_Dict.items()}
		Seq_as_str = str(Seq_ID_Seq_Dict[key])
		counter = forward_seqs_alignment_dict[key] - Seq_as_str.find(inverted_position_dict[position_of_match])
		Seq_as_list = list(Seq_as_str)
		# Initialize a counter to track where in a given sequence a given nucleotide occurs for output files
		seq_position = 0
		# Also creates a nested dictionary of first key as position relative to sequence, second key as sequence_ID, and value of position within sequence for a given sequence_ID
		for i in Seq_as_list:
			if counter in relative_position_of_nucleotide_in_seq_dict:
				if i in relative_position_of_nucleotide_in_seq_dict[counter]:
					relative_position_of_nucleotide_in_seq_dict[counter][i].update({key : seq_position}) 
					seq_position += 1
				else:
					relative_position_of_nucleotide_in_seq_dict[counter][i] = {key : seq_position}
			else:
				relative_position_of_nucleotide_in_seq_dict[counter] = {i : {key : seq_position}}
				seq_position += 1
			if counter in Contig_Dict:
				Contig_Dict[counter].append(key)
				Contig_Dict[counter].append(i)
				counter += 1
			else:
				Contig_Dict[counter] = [key, i]
				counter += 1
	return Contig_Dict, relative_position_of_nucleotide_in_seq_dict
