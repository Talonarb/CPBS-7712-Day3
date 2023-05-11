# Function to extract Query sequence as string and store the location of the sequence reads file for later use
import os.path
def Fasta_Reader(ReadDirectory, QueryDirectory, ReadFilename, QueryFilename):
    ReadFile = (os.path.join(ReadDirectory, ReadFilename))
    QueryFile = open(os.path.join(QueryDirectory, QueryFilename))
    FASTAquery = QueryFile.readlines()[1:]
    QueryFile.close
    # Strips any newline characters, converts query to string, and converts to a list excluding the first and last 2 characters, removing artifacts from reading/conversion
    FASTAquery = [item.strip() for item in FASTAquery]
    FASTAquery = str(FASTAquery)
    FASTAquery = list(FASTAquery[2:-2])
    return(ReadFile, FASTAquery)