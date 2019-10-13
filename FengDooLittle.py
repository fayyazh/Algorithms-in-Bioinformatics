import argparse
from NeedlemanWunsch import NeedlemanWunsch
import FileReader
from SumOfPairs import *
from Clustering import Clustering
import math

class FengDooLittle():

    #Variable Initialization
    def __init__(self, phylo_Tree, sequences, Scoring_Function, Gap_Penalty):
        self.phylo_Tree = phylo_Tree
        self.sequences = sequences
        self.Scoring_Function = Scoring_Function
        self.Gap_Penalty = Gap_Penalty
        self.errors = None
        self.scoringMatrix = []
        self.msa = []
        self.indices = []
        

    # Error check for correct letters in sequence and correct scoring matirx
    def errorCheck(self):
        errors = []

        #check all sequences
        nt = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
        for seq in self.sequences:
            for letter in seq:
                if (letter not in nt):
                    errors.append("Sequence contains invalid letter: " + letter)

        try:
            self.scoringMatrix = FileReader.ScoringMatrix(self, self.Scoring_Function)
        except FileNotFoundError:
            errors.append("file does not exist: " + self.Scoring_Function)
        except IOError:
            errors.append("file not readable: " + self.Scoring_Function)

        if (len(errors) == 0):
            self.errors = None
        else:
            self.errors = "\n".join(errors)

        return errors

    # Creates similarity matrix and then coversion to distance. Uses distance matrix for UPGMA/WPGMA and then eventually merge them to provide MSA.
    def MSA(self):

        sequences = self.sequences
        gap = self.Gap_Penalty
        Scoring_Function = self.Scoring_Function

        # build guide tree using UPGMA or WPGMA
        tree = Clustering(self.sequences, self.Gap_Penalty, self.Scoring_Function)
        #phyloTree(self, option)
        self.ChosenTree = tree.phyloTree(self.phylo_Tree)

        # Traversing the guided tree
        Clust_Nodes = tree.Sorted_Cluster(self.ChosenTree)

        msa = ["-" for i in range(len(sequences))]
        for nodes in Clust_Nodes:
            #In case of leaf node
            if nodes.distance == 0:
                msa[nodes.Leaf[0]] = sequences[nodes.Leaf[0]]
            elif len(nodes.Leaf) == 2:

                #aligning and merging sequences
                align = self.Pairwise_Alignment(sequences[nodes.Leaf[0]],
                                      sequences[nodes.Leaf[1]],
                                      Scoring_Function, gap)

                #replace gaps with *
                align.Aligned_Sequence[0] = align.Aligned_Sequence[0].replace('-', '*')
                align.Aligned_Sequence[1] = align.Aligned_Sequence[1].replace('-', '*')
                
                #update msa
                msa[nodes.Leaf[0]] = align.Aligned_Sequence[0]
                msa[nodes.Leaf[1]] = align.Aligned_Sequence[1]
            else:
                Group1 = nodes.child[0].Leaf
                Group2 = nodes.child[1].Leaf

                if len(Group1) > 1 and len(Group2) > 1:

                    #determine best pair and merging of groups
                    Align_Score = float('inf')
                    Best_Align = None
                    Best_Pair = (0, 0)
                    for i in range(len(Group1)):
                        for j in range(len(Group2)):
                            seq1 = msa[Group1[i]]
                            seq2 = msa[Group2[j]]
                            align = self.Pairwise_Alignment(seq1, seq2,
                                                  Scoring_Function, gap)
                            if align.finalAlignScore > Align_Score:
                                Align_Score = align.finalAlignScore
                                Best_Align = align.Aligned_Sequence
                                Best_Pair = (i, j)


                    # Deterimines indices of Gap in sequence
                    Seq1 = endresult[0]
                    Seq2 = endresult[1]
                    for i in range(len(Seq1)):
                        if Seq1[i] == "-":
                            self.indices.append(i)
                    Indice_Seq1 = self.indices

                    # merging 2 groups into 1 based on best pair
                    # alignment of other sequences in the group
                    for i in range(len(Group1)):
                        if i != Best_Pair[0]:
                            l = list(msa[Group1[i]])
                            for index in Indice_Seq1:
                                l.insert(index, '*')
                            msa[Group1[i]] = "".join(l)

                    for j in range(len(Seq2)):
                        if Seq1[j] == "-":
                            self.indices.append(j)
                    Indice_Seq2 = self.indices

                    for j in range(len(Group2)):
                        if j != Best_Pair[1]:
                            l = list(msa[Group2[j]])
                            for index in Indice_Seq2:
                                l.insert(index, '*')
                            msa[Group2[j]] = "".join(l)

                    # replace gaps with *
                    Best_Align[0] = Best_Align[0].replace('-', '*')
                    Best_Align[1] = Best_Align[1].replace('-', '*')
                    
                    #adjusting msa of the best pair
                    msa[Group1[Best_Pair[0]]] = Best_Align[0]
                    msa[Group2[Best_Pair[1]]] = Best_Align[1]
                else:
                    # merge sequence and group together
                    group = None
                    seq = None
                    if len(Group1) > 1:
                        group = Group1
                        seq = Group2
                    else:
                        group = Group2
                        seq = Group1

                    # determine best pair
                    Align_Score = float('-inf')
                    Best_Align = None
                    Best_Pair = (0, 0)
                    for i in range(len(group)):
                        align = self.Pairwise_Alignment(sequences[seq[0]], msa[group[i]],
                                              Scoring_Function, gap)
                        if align.finalAlignScore > Align_Score:
                            Align_Score = align.finalAlignScore
                            Best_Align = align.Aligned_Sequence
                            Best_Pair = (0, i)
    
                    #merge sequence into the group based on best pair
                    ind = Best_Align[1]
                    self.indices = []
                    for j in range(len(ind)):
                        if ind[j] == "-":
                            self.indices.append(j)
                    indices = self.indices

                    for i in range(len(group)):
                        if i != Best_Pair[1]:
                            l = list(msa[group[i]])
                            for index in indices:
                                l.insert(index, '*')
                            msa[group[i]] = "".join(l)
                            
                    # replace gaps with *
                    Best_Align[0] = Best_Align[0].replace('-', '*')
                    Best_Align[1] = Best_Align[1].replace('-', '*')
                    
                    #adjusting msa of the best pair
                    msa[seq[Best_Pair[0]]] = Best_Align[0]
                    msa[group[Best_Pair[1]]] = Best_Align[1]

        return msa
    

    # Determines pairwise alignment of sequences
    def Pairwise_Alignment(self, seq1, seq2, Scoring_Function, gap):
        align = NeedlemanWunsch()
        align.matrix = Scoring_Function
        align.Gap_Penalty = gap
        align.seq1 = seq1
        align.seq2 = seq2
        align.errorCheck()
        if align.error == "0":
            align.nw()
            return align
        else:
            print("Errors found:")
            Printer.PrintErrors(align.error)
            return None

    # Generates output
    def fdl_Results(self, MSA_Alignment):
        Result = []

        for s in MSA_Alignment:
            Result.append(self.output(s, 80))

        for i in range(len(Result[0])):
            for j in range(len(MSA_Alignment)):
                print("Seq " + str(j + 1) + " " + Result[j][i])
            print("\n")


    def output(self, string, length):
        subSection = math.ceil(len(string) / length)

        output = []
        for i in range(0, subSection):
            output.append(string[length * i: length * (i + 1)])

        return output

    def Errors(self, error):
        if (isinstance(error, list)):
            for err in error:
                print(err)
        elif (isinstance(error, str)):
            print(error)


def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastaFile", type=str)
    parser.add_argument("-g", "--Gap_Penalty", type=int, default=-8)
    parser.add_argument("-s", "--Scoring_Function", type=str, default="pam250", help="only file name. e.g. PAM250/BLOSUM62")
    parser.add_argument("-t", "--phylo_Tree", type=str, default="upgma")
    Pars = parser.parse_args()

    #reading fasta file and getting sequences
    sequences = []
    error = []
    sequences = FileReader.fastaFile(Pars.fastaFile)
    
    # feng doo little algorithm execution after error check
    FDL = FengDooLittle(Pars.phylo_Tree, sequences, Pars.Scoring_Function, Pars.Gap_Penalty)
    error = FDL.errorCheck()
    if len(error) == 0:
        MSA_Alignment = FDL.MSA()

        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("Feng-Doolittle Algorithm".center(80, ' ') + "\n")
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("INPUT SEQUENCES".center(80, ' ') + "\n")
        for i in range(len(sequences)):
            print("SEQUENCE " + str(i + 1) + ": " + "\n" + sequences[i])
            print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("GAP PENALTY".center(80, ' '))
        print("Gap Penalty: " + str(FDL.Gap_Penalty))
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("SCORING FUNCTION".center(80, ' '))
        print("Scoring Function: " + FDL.Scoring_Function)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("ALIGNMENT".center(80, ' ') + "\n")
        FDL.fdl_Results(MSA_Alignment)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        SOP = Sum_Of_Pairs([], FDL.Scoring_Function)
        print("SUM OF PAIR SCORE".center(80, ' '))
        print("Sum of Pair Score: " + str(SOP.ComputeScore(MSA_Alignment)))
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    else:
        print("Errors found:")
        FDL.Errors(error)

if __name__ == "__main__":
    Main()
