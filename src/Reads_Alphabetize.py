# Originally intended to sort sequences alphabetically, currently just used to create a list of sequences to pass to another function
from . import Summary_Info
def Reads_Alphabetize(file, Fasta_query):
    import numpy as np
    import pandas as pd
    reads_df = pd.read_csv(file, sep = "|", names = ['data'])
    reads_df = pd.DataFrame({'Seq_ID': reads_df.data.iloc[::2].to_numpy(), 'Seq_Read': reads_df.data.iloc[1::2].to_numpy()}) 
    #print(reads_df.head())
    #reads_df = reads_df.sort_values(by=['Seq_Read'])
    #print(reads_df.head())
    #Seq_ID_list = reads_df.Seq_ID.values.tolist()
    Seq_read_list = reads_df.Seq_Read.values.tolist()
    return Summary_Info.Summary_Info(Seq_read_list, Fasta_query)
    #Merged_list = Merge_for_file_generation(Seq_ID_list, Seq_read_list)
    #sorted_fasta = open("SORTED.fasta", "w")
    #for i in Merged_list:
    #    sorted_fasta.write(i)
    #    sorted_fasta.write('\n')
    #sorted_fasta.close