import sys
import os
import unittest
from src import Seq_Comp, Fasta_Filter, Reads_Alphabetize, Complement_Conversion, Position_Dict_Maker, Alignment_of_Seqs, Contig_Dict_Maker, Contig_Assembler

class Test_Seq_Comp(unittest.TestCase):
    def test_comp(self):
        #Tests if the composition of a given item is deconstructed into a set of overlapping 3-character units
        test_query = "ABC123*/+"
        self.assertEqual(Seq_Comp.Seq_Comp(test_query, 3), {('A','B','C'), ('B','C','1'), ('C','1','2'), ('1','2','3'), ('2','3','*'), ('3','*','/'),('*','/','+')})

class Test_Fasta_Filter(unittest.TestCase):
    def test_Fasta_Filter(self):
        fasta_to_check = os.path.join(os.path.dirname(__file__), 'TEST.fasta')
        # Original test query composition "THISISATESTFOR", fasta file contains reads that correspond except for one ">4:2:3", which would ony overlap for 3 characters
        query_for_ff = {('T','H','I','S'),('H','I','S','I'),('I','S','I','S'),('S','I','S','A'),('I','S','A','T'),('S','A','T','E'),('A','T','E','S'),('T','E','S','T'),('E','S','T','F'),('S','T','F','O'),('T','F','O','R')}
        sequences_to_keep = Fasta_Filter.Fasta_filter(fasta_to_check, query_for_ff,4)
        self.assertEqual(sequences_to_keep, [">1:2:3",">2:2:3",">3:2:3"])

class Test_Reads_Alphabetize(unittest.TestCase):
    def test_Reads_Alphabetize(self):
        
        return

if __name__ == '__main__':
    unittest.main()