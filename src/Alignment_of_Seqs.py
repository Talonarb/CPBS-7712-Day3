def Alignment_of_Seqs(Query_dict, seqs_to_keep, ReadFile, suggested_bp_overlap):
    read_alignments_to_query = {}
    Seq_ID_Seq_Dict = {}
    Seq_25mer_aligned = {}
    with open(ReadFile) as f:
        line = next(f)
        while line:
            try:
                if line[0] == ">":
                    Seq_ID = line
                    Stripped_Seq_ID = Seq_ID.strip()
                    line=next(f)
                    if Stripped_Seq_ID not in seqs_to_keep:
                        line=next(f)
                    else:
                        Seq_ID_Seq_Dict[Stripped_Seq_ID] = line.strip()
                        read = list(line)
                        read_part = []
                        for i in range(len(read)):
                            read_part.append(read[i])
                            i+=1
                            if len(read_part) == suggested_bp_overlap:
                                str_of_olap = ''.join(map(str,read_part))
                                if str_of_olap in Query_dict.keys():
                                    read_alignments_to_query[Stripped_Seq_ID] = Query_dict[str_of_olap]
                                    Seq_25mer_aligned[Stripped_Seq_ID] = str_of_olap
                                read_part = read_part[1:suggested_bp_overlap]
                        line = next(f)
            except StopIteration as e:
                break
    return read_alignments_to_query, Seq_ID_Seq_Dict, Seq_25mer_aligned