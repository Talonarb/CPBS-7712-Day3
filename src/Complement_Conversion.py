# Simple helper function to convert DNA (as list, with elements of list as nucleotides)
# to the reverse complement of the strand
# Converts sequence to its reverse complement when sequence is given as a list of nucleotides 
def complement_conversion(seq_to_convert):
    converted_seq = []
    for i in seq_to_convert:
        if i == "A":
            converted_seq.append("T")
        elif i == "T":
            converted_seq.append("A")
        elif i == "C":
            converted_seq.append("G")
        elif i == "G":
            converted_seq.append("C")
        elif i != "A" or "C" or "T" or "G":
            print("The sequence does not consist of only nucleotides, please double check your input")
            break
    return(converted_seq)

