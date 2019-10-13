import argparse
from NeedlemanWunsch import NeedlemanWunsch
import FileReader
from SumOfPairs import *
import math
import copy

class Phylo_Tree:
    def __init__(self, child, distance = 0, parentDist = 0):
        self.Leaf = []
        self.child = child
        self.distance_Parent = parentDist
        self.distance = distance

class Clustering:
    def __init__(self, sequences, gap, scoring, root = None):
        self.sequences = sequences
        self.Tree = None
        self.Gap_Penalty = gap
        self.Scoring_Function = scoring
        self.scoringMatrix = []
        self.root = []
        self.errors = []

        if(root is None):
            self.root = []
            for i in range (1, len(sequences) + 1):
                self.root.append("Seq" + str(i))
        else:
            self.root = root

    #Error Check
    def errorCheck(self):
        err = []

        #check possible errors
        nt = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
        for seq in self.sequences:
            for letter in seq:
                if (letter not in nt):
                    err.append("Sequence contains invalid letter: " + letter)

        try:
            self.scoringMatrix = FileReader.ScoringMatrix(self, self.Scoring_Function)
        except FileNotFoundError:
            err.append("file is not found: " + self.Scoring_Function)
        except IOError:
            err.append("file cannot be opened: " + self.Scoring_Function)

        if (len(err) == 0):
            self.errors = None
        else:
            self.errors = "\n".join(err)

        return err

    def phyloTree(self, option):

        #compute similarity matrix
        similarityMatrix = self.Similarity_Matrix()

        #convert similarity to distance
        distanceMatrix = self.Similarity_Distance(similarityMatrix)

        clusters = self.PhyloCluster(distanceMatrix, option)

        self.Tree = clusters[0]
        return clusters[0]

    def PhyloCluster(self, distanceMatrix, option):
        clusters = []
        root = self.root
        for i in range(len(root)):
            c = Phylo_Tree([root[i]], 0)
            c.Leaf.append(i)
            clusters.append(c)

        distMat = copy.deepcopy(distanceMatrix)

        # combining two nearest cluster with minimum distance
        while (len(clusters) > 1):
            minPair = self.Min_Distance(distMat)

            #calculate distance to parent node
            clusters[minPair[0]].distance_Parent = distMat[minPair[0]][minPair[1]] / 2 - clusters[minPair[0]].distance
            clusters[minPair[1]].distance_Parent = distMat[minPair[0]][minPair[1]] / 2 - clusters[minPair[1]].distance

            #join the minimum distance pair
            c = Phylo_Tree([clusters[minPair[0]], clusters[minPair[1]]],
                        distMat[minPair[0]][minPair[1]] / 2)
            c.Leaf.extend(clusters[minPair[0]].Leaf)
            c.Leaf.extend(clusters[minPair[1]].Leaf)
            clusters.append(c)

            if (option == "upgma"):
                # delete in descending index order
                for index in sorted(minPair, reverse=True):
                    # delete from clusters list
                    del clusters[index]
                    # delete row from distMat
                    del distMat[index]
                    # delete col from distMat
                    for row in distMat:
                        del row[index]
                # update distance matrix
                distMat = self.Update_UPGMA(clusters, distMat, distanceMatrix)

            elif (option == "wpgma"):
                distMat = self.Update_WPGMA(clusters, distMat, minPair)

                # delete in descending index order
                for index in sorted(minPair, reverse=True):
                    # delete from clusters list
                    del clusters[index]
                    # delete row from distMat
                    del distMat[index]
                    # delete col from distMat
                    for row in distMat:
                        del row[index]

        return clusters

    # Create Similarity matrix using needlemanwunsch pairwise alignment to determine distance matrix
    def Similarity_Matrix(self):
        sequences = self.sequences
        gap = self.Gap_Penalty
        Scoring_Function = self.Scoring_Function
        sop = Sum_Of_Pairs(sequences, self.Scoring_Function)

        similarityMatrix = [[None for i in range(len(sequences))]for j in range(len(sequences))]
        
        for i in range(len(sequences)):
            for j in range(i, len(sequences)):
                if i != j:
                    nw = self.Pairwise_Alignment(sequences[i], sequences[j],
                                          Scoring_Function, gap)
                        
                    similarityMatrix[i][j] = nw.finalAlignScore
                    similarityMatrix[j][i] = nw.finalAlignScore
                else:
                    similarityMatrix[i][i] = sop.ComputeScore([sequences[i],
                                                                 sequences[i]])

        return similarityMatrix


    # Converts Simalirty Matrix into Distance to build tree
    def Similarity_Distance(self, similarityMatrix):
        gap = self.Gap_Penalty
        Scoring_Function = self.Scoring_Function
        sequences = self.sequences

        distanceMatrix = [[None for i in range(len(sequences))]for j in range(len(sequences))]
        
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                seq1 = sequences[i]
                seq2 = sequences[j]

                nw = self.Pairwise_Alignment(seq1, seq2, Scoring_Function, gap)
                sRand = self.S_Rand(nw)

                #S_Max = self.S_Max(nw)
                S_Max = (similarityMatrix[i][i] + similarityMatrix[j][j]) / 2
                S_Eff = (similarityMatrix[i][j] - sRand) / (S_Max - sRand)

                distance = - (math.log(S_Eff))

                distanceMatrix[i][j] = distance
                distanceMatrix[j][i] = distance

            distanceMatrix[i][i] = 0

        return distanceMatrix



    # Update distance matrix in case of UPGMA
    def Update_UPGMA(self, clusters, distMat, Initial_Matrix):
        distanceMatrix = Initial_Matrix

        distMat.append([None for i in range(len(clusters))])
        
        #last cluster is always the new cluster (because list.append)
        lastIndex = len(clusters) - 1
        newCluster = clusters[lastIndex]

        #for each cluster except new cluster
        for i in range(len(clusters[:-1])):
            cluster = clusters[i]
            numerator = 0
            denominator = 0

            #calculate distance to new cluster
            for Old_Leaf in cluster.Leaf:
                for New_Leaf in newCluster.Leaf:
                    #get distances to each other                    
                    numerator += distanceMatrix[Old_Leaf][New_Leaf]
                    denominator += 1
                
            #take arithmatic mean
            mean = numerator / denominator
            distMat[lastIndex][i] = mean
            distMat[i].append(mean)

        #distance of new cluster to itself is 0
        distMat[lastIndex][lastIndex] = 0

        return distMat


    # Update matrix in case of WPGMA
    def Update_WPGMA(self, clusters, distMat, Best_Pair):
        distMat.append([None for i in range(len(clusters))])
        lastIndex = len(clusters) - 1

        #for each cluster except new cluster
        for i in range(len(clusters[:-1])):
            if(i != Best_Pair[0]) and (i != Best_Pair[1]):
                cluster = clusters[i]
                
                #calculate distance to new cluster
                dist1 = distMat[Best_Pair[0]][i]
                dist2 = distMat[Best_Pair[1]][i]
                newDist = (dist1 + dist2) / 2

                distMat[i].append(newDist)
                distMat[lastIndex][i] = newDist

        distMat[lastIndex][lastIndex] = 0

        return distMat

    # Search Minimum distance in the distance matrix
    def Min_Distance(self, matrix):
        minVal = float('inf')
        minPair = (0, 0)

        for row in range(len(matrix)):
            for col in range(len(matrix)):
                if (matrix[row][col] < minVal and matrix[row][col] > 0):
                    minVal = matrix[row][col]
                    minPair = (row, col)

        return minPair

    # Sorted Clusers
    def Sorted_Cluster(self, Tree, Decr_Dist = False):
        nodes = [Tree]
        for C in Tree.child:
            if (type(C) is Phylo_Tree):
                nodes.extend(self.Sorted_Cluster(C))
        
        return sorted(nodes, key=lambda c: c.distance, reverse = Decr_Dist)

    # Calculates S_Rand
    def S_Rand(self, nw):
        alignedSeqA = nw.Aligned_Sequence[0]
        alignedSeqB = nw.Aligned_Sequence[1]
        gapcost = self.Gap_Penalty
        noOfGaps = 0

        # count no of gaps
        for i in range(len(alignedSeqA)):
            if alignedSeqA[i] == '-':
                noOfGaps += 1
        for i in range(len(alignedSeqB)):
            if alignedSeqB[i] == '-':
                noOfGaps += 1

        # build frequency dict
        dictA = self.Char_Frequency(alignedSeqA)
        dictB = self.Char_Frequency(alignedSeqB)

        summation = 0
        for keyA in dictA:
            for keyB in dictB:
                if keyA != '-' and keyB != '-':
                    score = int(nw.scoringMatrix[keyA][keyB])
                    freqA = dictA[keyA]
                    freqB = dictB[keyB]
                    summation += score * freqA * freqB

        sRand = (1 / len(alignedSeqA)) * summation - noOfGaps
        # sRand = (summation - (noOfGaps * gapcost)) / len(alignedSeqA)

        return sRand

    def S_Max(self, nw):
        alignedSeqA = nw.variables.Aligned_Sequence[0]
        alignedSeqB = nw.variables.Aligned_Sequence[1]

        dictA = self.Char_Frequency(alignedSeqA)
        dictB = self.Char_Frequency(alignedSeqB)

        S_Max_aa = 0
        for keyA in dictA:
            if keyA != '-':
                score = nw.GetScore(keyA, keyA)
                S_Max_aa += score

        S_Max_bb = 0
        for keyB in dictB:
            if keyB != '-':
                score = nw.GetScore(keyB, keyB)
                S_Max_bb += score

        S_Max = (S_Max_aa + S_Max_bb) / 2

        return S_Max

    # Frequency of Characters
    def Char_Frequency(self, sequence):
        frequency = {}
        for char in sequence:
            if char not in frequency:
                frequency[char] = 1
            else:
                frequency[char] = frequency[char] + 1

        return frequency


    # pairwise alignment for similarity matrix
    def Pairwise_Alignment(self, seq1, seq2, Scoring_Function, gap):
        nw = NeedlemanWunsch()
        nw.matrix = Scoring_Function
        nw.Gap_Penalty = gap
        nw.seq1 = seq1
        nw.seq2 = seq2
        nw.errorCheck()
        if nw.error == "0":
            nw.nw()
            return nw
        else:
            return None

    def Newick_Format(self, Tree = None, root = None):
        if Tree is None:
            Tree = self.Tree
        if root is None:
            root = self.root
            
        Newick = {}
        cluster = copy.deepcopy(Tree)
        rootDistance = cluster.distance

        Clust_Nodes = self.Sorted_Cluster(cluster)

        for nodes in Clust_Nodes:
            if len(nodes.Leaf) == 1:
                Newick[nodes] = root[nodes.Leaf[0]] + ":" + str(nodes.distance_Parent)
            else:
                Newick[nodes] = "(" + Newick[nodes.child[0]] + ", " + Newick[nodes.child[1]] \
                + "):" + str(nodes.distance_Parent)

        print(Newick[cluster])

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
    parser.add_argument("-tree", "--upgma", type=str, default="upgma")
    parser.add_argument("fastaFile", type=str)
    parser.add_argument("-g", "--Gap_Penalty", type=int, default=-8)
    parser.add_argument("-s", "--Scoring_Function", type=str, default="pam250", help="File extension not required, only name. e.g. PAM250/Blosums62")

    # reading Input from a fasta file
    Pars = parser.parse_args()
    sequences = []
    sequences = FileReader.fastaFile(Pars.fastaFile)
    
    # Perform clustering according to input after error check
    errors = []
    cluster = Clustering(sequences, Pars.Gap_Penalty, Pars.Scoring_Function)
    errors = cluster.errorCheck()
    if len(errors) == 0:
        option = cluster.phyloTree(Pars.upgma)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("CLUSTERING (UPGMA/WPGMA)".center(80, ' ') + "\n")
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("SEQUENCES".center(80, ' ') + "\n")
        for i in range(len(sequences)):
            print("SEQUENCE " + str(i + 1) + ": " + "\n" + sequences[i])
            print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("GAP PENALTY".center(80, ' '))
        print("Gap Penalty: " + str(cluster.Gap_Penalty))
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("SCORING FUNCTION".center(80, ' '))
        print("Scoring Function: " + cluster.Scoring_Function)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("Result Newick Format:\n".center(80, ' ') + "\n")
        cluster.Newick_Format(option)
        #cluster.Newick_Format(Chosen_Tree)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")

    else:
        print("Errors found:")
        cluster.Errors(errors)

if __name__ == "__main__":
    Main()
