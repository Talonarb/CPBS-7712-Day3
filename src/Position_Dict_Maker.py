def position_dict_maker(Query, bp_overlap_minimum):
    composition = []
    position_dict = {}
    position_counter = 0
    for i in range(len(Query)):
        composition.append(Query[i])
        i+=1
        if len(composition) == bp_overlap_minimum:
            str_of_olap = ''.join(map(str,composition))
            position_dict[str_of_olap] = position_counter
            position_counter +=1
            composition = composition[1:bp_overlap_minimum]
            i+=1
    return(position_dict)
