# Query Sequence Assisted Assembly

## Goal
Given a Query in FASTA format as well as reads in a FASTA format, assemble the longest contig containing the Query.

## Description
This program takes two FASTA files as input, a Query FASTA and a sequences FASTA and outputs a contig assembled from the sequences containing the Query. 

Some calculations take place to choose a minimum overlap for assembly, and sequence reads are filtered out until those remaining overlap by a minimum of 
the calculated number. Following this, the positions of the overlap, the sequence ID and it's respective sequence, and the composition of the overlap are
used to assemble an initial contig based on consensus of sequences. Folowing this steps, any gaps or non-contiguous portions are filtered out, and the final
contig is assembled. Finally, the alignment tsv is output showing the overlaps of all sequences to the query, although it should be noted coverage includes 
the ability for any overlap at position 25 to be a mismatch. The program can take user inputs, but defaults to the a file titled FASTA where it expects to find
QUERY.fasta and READS.fasta. 

 

## Install
- pandas~=1.4.1

## Usage

### Python
```python 
# Author: Talon Arbuckle

ReadDirectory = (input("Enter the directory of the FASTA reads file: ")) or (os.path.join(sys.path[0])+"\FASTA")
QueryDirectory = (input("Enter the directory of the QUERY query file: ")) or (os.path.join(sys.path[0])+"\FASTA")
ReadFilename = (input("Enter the filename of the FASTA reads file: ")) or ("READS.fasta")
QueryFilename = (input("Enter the filename of the  reads file: ")) or ("QUERY.fasta") 
def main():
    # Read in files
    Read_file, Fasta_query = Fasta_Reader.Fasta_Reader(ReadDirectory, QueryDirectory, ReadFilename, QueryFilename)
    # Calculate the overlap suggested based on query length and average sequence length
    suggested_bp_overlap_minimum = Reads_Alphabetize.Reads_Alphabetize(Read_file, Fasta_query)
    # Break the query down to into overlapping parts of the suggesteed overlap, tracking positions of each
    Query_comp = (Seq_Comp.Seq_Comp(Fasta_query, suggested_bp_overlap_minimum))
    Position_Dict = Position_Dict_Maker.position_dict_maker(Fasta_query, suggested_bp_overlap_minimum)
    # Filter the sequences using the list of diction
    Sequences_to_keep = Fasta_Filter.Fasta_filter(Read_file, Query_comp, suggested_bp_overlap_minimum)
    # Create several dictionaries for use in next steps, including which position the sequences align to, the sequence ID and sequence itself, and 
    # the composition of which broken down section of the query the sequence matches
    Forward_Seqs_Alignment_Dict, Seq_ID_Seq_Dict, Seq_suggeseted_mer_aligned = Alignment_of_Seqs.Alignment_of_Seqs(Position_Dict, Sequences_to_keep, Read_file, suggested_bp_overlap_minimum)
    # Code to do so for the complement of the query, no results, so no next steps taken
    Reverse_Comp_of_Query = Complement_Conversion.complement_conversion(Fasta_query)
    Query_Comp_Reverse = (Seq_Comp.Seq_Comp(Reverse_Comp_of_Query, suggested_bp_overlap_minimum))
    Sequences_to_keep_reverse = Fasta_Filter.Fasta_filter(Read_file, Query_Comp_Reverse, suggested_bp_overlap_minimum)
    # Creation of a contig dict with gaps as well as the sequence of origin for every base used in the first contig
    Contig_Dict, Relative_Position_of_Nucleotide_in_Sequence_dict = Contig_Dict_Maker.Contig_Dict_Maker(Forward_Seqs_Alignment_Dict, Seq_ID_Seq_Dict, Position_Dict)
    # Filters out gaps from first assembly of overlap, assembles final contig
    Final_Seq_Comp_Dict, Final_contig = Contig_Assembler.contig_assembly(Contig_Dict, Relative_Position_of_Nucleotide_in_Sequence_dict)
    # Converts the query to a string for printing, and finds the position of the query in the final assembles contig
    Fasta_String = ''.join(map(str, Fasta_query))
    query_position_in_contig = (Final_contig.find(Fasta_String))
    # Outfile for alignment and contig
    tsv_generator.comp_tsv_generator(Final_Seq_Comp_Dict, Seq_suggeseted_mer_aligned, Seq_ID_Seq_Dict, Position_Dict, Fasta_String, suggested_bp_overlap_minimum)

main()
```
### Command Line
$ python main.py 
Enter the directory of the FASTA reads file: 
Enter the directory of the QUERY query file: 
Enter the filename of the FASTA reads file: 
Enter the filename of the  reads file: 
```
$ python main.py 


## Input
- Query and Sequence reads files in FASTA format

## Output
- tsv file with alignment information of final contig
- final contig as FASTA file
