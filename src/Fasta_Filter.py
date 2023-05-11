# Function to read individual sequences from a fasta file, which outputs sequence IDs that have a sequence in which the set K-mers that comprise it have at least one interesection
# with the set of K-mers generated from the query sequence.  
from . import Seq_Comp
def Fasta_filter(ReadFile, Query_comp, bp_overlap_minimum):
    Sequences_to_keep = []
    with open(ReadFile) as f:
        line = next(f)
        while line:
            try:
                if line[0] == ">":
                    Seq_ID = line
                    line = next(f)
                else:
                    Sequence_comp = Seq_Comp.Seq_Comp(line, bp_overlap_minimum)
                    if Query_comp.intersection(Sequence_comp):
                        Stripped_Seq_ID = Seq_ID.strip()
                        Sequences_to_keep.append(Stripped_Seq_ID)
                    line=next(f)
            except StopIteration as e:
                break
    return(Sequences_to_keep) 