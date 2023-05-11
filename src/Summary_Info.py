def Summary_Info(Read_list, QUERY):
    from functools import reduce
    from math import floor
    avg_read = reduce(lambda x, y: x + y, map(len, Read_list))/len(Read_list)
    #print("The average read length is", avg_read)
    #print("The length of the query is", len(QUERY))
    #print("The length of query/average length is", (len(QUERY)/avg_read))
    suggested_bp_overlap_minimum = floor((avg_read)/((len(QUERY)/avg_read)))
    return suggested_bp_overlap_minimum