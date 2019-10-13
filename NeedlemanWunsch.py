import FileReader
import argparse
import random
import math

class NeedlemanWunsch:
    def __init__(self):
        self.scoringMatrix = None
        self.finalAlignScore = None
        self.seq1 = None
        self.seq2 = None
        self.Gap_Penalty = - 10
        self.Scoring_Matrix = "pam250"
        self.Filled_Matrix = None
        self.Aligned_Sequence = None
        self.error = "0"
        self.finalAlignScore = None
        self.scoringMatrix = None
        self.value = 0
        self.source = []

    def __repr__(self):
        return str(self.source)

    # Verifies the input and errors list is updated
    def errorCheck(self):
        errors = []
        
        if not isinstance(self.Gap_Penalty, int):
            errors.append("Gap cost value should be integer")

        nt = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
        for seq in self.seq1:
            for letter in seq:
                if (letter not in nt):
                    errors.append("Sequence contains invalid letter: " + letter)

        nt = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
        for seq in self.seq2:
            for letter in seq:
                if (letter not in nt):
                    errors.append("Sequence contains invalid letter: " + letter)

        try:
            self.scoringMatrix = FileReader.ScoringMatrix(self, self.Scoring_Matrix)
        except FileNotFoundError:
            errors.append("file is not found: " + self.Scoring_Matrix)
        except IOError:
            errors.append("file cannot be opened: " + self.Scoring_Matrix)

        if (len(errors) == 0):
            self.errors = None
        else:
            self.errors = "\n".join(errors)

        return errors

     # Core work of needlemanwunsch is being executed here
    def nw(self):
        self.scoringMatrix = FileReader.ScoringMatrix(self, self.Scoring_Matrix)
        self.Filled_Matrix= self.Matrix_Filling()
        self.Aligned_Sequence = self.Traceback()

        return self.Aligned_Sequence

    # Matrix is being initialized for later use
    def Matrix_Initialization(self):

        Initialized_matrix = [[NeedlemanWunsch() for i in range(len(self.seq1) + 1)] for i in range(len(self.seq2) + 1)]

        #init horizontal
        hor = 1
        for alphabet in self.seq1:
            Initialized_matrix[0][hor].value = Initialized_matrix[0][hor - 1].value + self.Gap_Penalty
            Initialized_matrix[0][hor].source.append(str(0) + ',' + str(hor - 1) + ',LEFT')
            hor += 1

        #init vertical
        ver = 1
        for alphabet in self.seq2:
            Initialized_matrix[ver][0].value = Initialized_matrix[ver - 1][0].value + self.Gap_Penalty
            Initialized_matrix[ver][0].source.append(str(ver - 1) + ',' + str(0) + ',UP')
            ver += 1
            
        return Initialized_matrix

    # Filling the matrix and final alignment score is saved
    def Matrix_Filling(self):

        matrix = self.Matrix_Initialization()
        
        hor_cell = 1
        ver_cell = 1

        for B in self.seq2:
            for A in self.seq1:
                Left = matrix[ver_cell][hor_cell - 1].value + self.Gap_Penalty
                Up = matrix[ver_cell - 1][hor_cell].value + self.Gap_Penalty
                Diagonal = matrix[ver_cell - 1][hor_cell - 1].value + int(self.scoringMatrix[A][B])

                value = max(Up, Left, Diagonal)
                matrix[ver_cell][hor_cell].value = value

                if (value == Left):
                    matrix[ver_cell][hor_cell].source.append(str(ver_cell) + ',' + str(hor_cell - 1) + ',LEFT')
                if(value == Up):
                    matrix[ver_cell][hor_cell].source.append(str(ver_cell - 1) + ',' + str(hor_cell) + ',UP')
                if(value == Diagonal):
                    matrix[ver_cell][hor_cell].source.append(str(ver_cell - 1) + ',' + str(hor_cell - 1) + ',DIAGONAL')

                hor_cell += 1
            ver_cell += 1
            hor_cell = 1

        # final alignment score
        self.finalAlignScore = matrix[len(self.seq2)][len(self.seq1)].value

        return matrix

    # value from up and left requires insertion of gap and incase of diagonal match/mismatch
    def Traceback(self):

        Column = len(self.seq1)
        Row = len(self.seq2)
        SequenceA = Column - 1
        SequenceB = Row - 1
        Const_SeqA = []
        Const_SeqB = []

        while Row != 0 or Column != 0:
            # Determine source of the cell value
            source = self.Filled_Matrix[Row][Column].source
            comp = source[random.randint(0, len(source) - 1)].split(',')
            Row = int(comp[0])
            Column = int(comp[1])

            if comp[2] == "UP":
                Const_SeqA.insert(0, '-')
                Const_SeqB.insert(0, self.seq2[SequenceB])
                SequenceB -= 1
            elif comp[2] == "LEFT":
                Const_SeqA.insert(0, self.seq1[SequenceA])
                Const_SeqB.insert(0, '-')
                SequenceA -= 1
            else:
                Const_SeqA.insert(0, self.seq1[SequenceA])
                Const_SeqB.insert(0, self.seq2[SequenceB])
                SequenceA -= 1
                SequenceB -= 1

        alignedSeqA = ''.join(Const_SeqA)
        alignedSeqB = ''.join(Const_SeqB)

        return [alignedSeqA, alignedSeqB]


    def NWAlign(self, Result):
        if (not isinstance(Result[0], list)):
            resultA = self.output(Result[0], 80)
            resultB = self.output(Result[1], 80)
            seperator = self.Seperator(Result[0], Result[1])
            resultFormat = self.output(seperator, 80)

            for i in range(0, len(resultA)):
                a = i + 1
                print("Alignment : ", a, "\n")
                print(resultA[i])
                print(resultFormat[i])
                print(resultB[i])
                print("...............................................................................")
                print("\n")
        else:
            i = 1
            for result in Result:
                print("Result " + str(i))
                NWAlign(result)
                print("\n")
                i += 1
    def Seperator(self, seqA, seqB):
        Sep = []
        i = 0
        for letter in seqA:
            if (letter == '-' or seqB[i] == '-'):
                Sep.append(' ')
            elif (letter == seqB[i]):
                Sep.append('*')
            elif (letter != seqB[i]):
                Sep.append(':')
            i += 1

        return "".join(Sep)

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
    parser.add_argument("-g", "--gapCost", type=int, default=-4)
    parser.add_argument("-s", "--Scoring_Function", type=str, default="pam250", help="only file name. e.g. PAM250/BLOSUM62")

    Pars = parser.parse_args()

    #reading fasta file and getting sequences
    seq1 = []
    seq2 = []
    error = []
    sequences = FileReader.fastaFile(Pars.fastaFile)
    #if more than 2 sequence is found, take the first 2
    seq1 = sequences[0]
    seq2 = sequences[1]


    #computing alignment
    nw = NeedlemanWunsch()
    nw.seq1 = seq1
    nw.seq2 = seq2
    nw.Gap_Penalty = Pars.gapCost
    nw.matrix = Pars.Scoring_Function
    nw.errorCheck()
    error = nw.error
    if(error == "0"):
        alignmentResult = nw.nw()


        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("NEEDLEMANWUNSCH ALGORITHM".center(80, ' ') + "\n")
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("INPUT SEQUENCES".center(80, ' ') + "\n")
        print("Sequence 1: " + seq1)
        print("Sequence 2: " + seq2)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("GAP PENALTY".center(80, ' '))
        print("Gap Cost: " + str(nw.Gap_Penalty))
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("SCORING FUNCTION".center(80, ' '))
        print("Scoring Matrix: " + nw.matrix)
        print("--------------------------------------------------------------------------------")
        print("Alignment Result:".center(80, ' ') + "\n")
        nw.NWAlign(alignmentResult)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("Alignment Score:".center(80, ' ') + "\n")
        print("Score: " + str(nw.finalAlignScore))
        print("--------------------------------------------------------------------------------")
    else:
        print("Errors found:")
        nw.Errors(error)


if __name__ == "__main__":
    Main()

