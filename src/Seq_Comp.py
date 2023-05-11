# This helper function iterates over a sequence, creating a set of all unique K-mers based on a suggested overlap calculated with average read length and query length. 
def Seq_Comp(Query, bp_overlap_minimum):
    Query_list = [x for x in Query]
    composition = []
    Composition_list = []
    for i in range(len(Query_list)):
        composition.append(Query_list[i])
        i+=1
        if len(composition) == bp_overlap_minimum:
            Composition_list.append(composition)
            composition = composition[1:bp_overlap_minimum]
            i+=1
    # This converts the list of lists to a set, which is faster to compare to than a list of lists
    Composition_set = {tuple(x) for x in Composition_list}
    return(Composition_set)    