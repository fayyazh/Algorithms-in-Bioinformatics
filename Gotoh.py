import argparse
import random
import FileReader
import sys
import math

class matrixConst:
    def __init__(self, val = 0):
        self.value = val
        self.source = []

    def __repr__(self):
        return str(self.source)

class Gotoh:
    def __init__(self):
        self.errors = None
        self.sequenceA = None
        self.sequenceB = None
        self.gapPenality = -4
        self.gapExtended = -1
        self.IP_Matrix = "pam250"
        self.filledMatrices = None
        self.alignedSequence = None
        self.finalAlignScore = None
        self.scoringMatrix = None



    # verifies the inputs and Errors list is updated.
    def errorCheck(self):
        errors = []

        if not isinstance(self.gapPenality, int):
            errors.append("Gap cost value should be integer")
        if not isinstance(self.gapExtended, int):
            errors.append("Gap extension value should be integer")
                
        try:
            self.scoringMatrix = FileReader.ScoringMatrix(self, self.IP_Matrix)
        except FileNotFoundError:
            errors.append("file is not found: " + filename)
        except IOError:
            errors.append("file cannot be opened: " + filename)

        except ValueError as Err:
            errors.append("values are not integer: " + Err)

        except:
            err = sys.exc_info()[0]
            errors.append("Error: " + str(err))

        if(len(errors) == 0):
            self.errors = "0"
        else:
            self.errors = "\n".join(errors)

    # core work of needlemanWunsch alignment is being executed here.
    def gotoh(self):
        self.scoringMatrix = FileReader.ScoringMatrix(self, self.IP_Matrix)
        self.filledMatrices = self.matrixFilling()
        self.alignedSequence = self.Traceback()
        return self.alignedSequence

    # The matrix is being initialized for later use to fill the matrix as per NeedlemanWunsch Algorithm
    def matrixInitialization(self):
        matrix_P = [[matrixConst() for i in range(len(self.sequenceA) + 1)] for i in range(len(self.sequenceB) + 1)]
        matrix_Q = [[matrixConst() for i in range(len(self.sequenceA) + 1)] for i in range(len(self.sequenceB) + 1)]
        matrix_D = [[matrixConst() for i in range(len(self.sequenceA) + 1)] for i in range(len(self.sequenceB) + 1)]

        for cVal in matrix_P[0][1:]:
            cVal.value = float('-inf')

        for rVal in matrix_P[1:]:
            rVal[0].value = None
            
        for cVal in matrix_Q[0][1:]:
            cVal.value = None

        for rVal in matrix_Q[1:]:
            rVal[0].value = float('-inf')

        i = 1
        for cVal in matrix_D[0][1:]:
            cVal.value = self.afineGap(i)
            cVal.source.append("2,0," + str(i - 1) + ",LEFT")
            i += 1
            
        i = 1
        for rVal in matrix_D[1:]:
            rVal[0].value = self.afineGap(i)
            rVal[0].source.append("2," + str(i - 1) + ",0,UP")
            i += 1

        return [matrix_P, matrix_Q, matrix_D]

    # Afine Gap calculation
    def afineGap(self, Gaps):
        return self.gapPenality + Gaps * self.gapExtended

    # Fill the matrix and the final alignment score is saved.
    def matrixFilling(self):
        matrices = self.matrixInitialization()
        verCell = 1
        horCell = 1

        for rVal in matrices[0][1:]:
            for cVal in matrices[0][verCell][1:]:

                left_Q = matrices[1][verCell][horCell - 1].value + self.gapExtended
                left_D = matrices[2][verCell][horCell - 1].value + self.afineGap(1)
                value_Q = max(left_Q, left_D)
                matrices[1][verCell][horCell].value = value_Q

                if(value_Q == left_Q):
                    matrices[1][verCell][horCell].source.append("1," + str(verCell) + "," + str(horCell - 1) + ",LEFT")
                if(value_Q == left_D):
                    matrices[1][verCell][horCell].source.append("2," + str(verCell) + "," + str(horCell - 1) + ",LEFT")

                up_P = matrices[0][verCell - 1][horCell].value + self.gapExtended
                up_D = matrices[2][verCell - 1][horCell].value + self.afineGap(1)
                value_P = max(up_P, up_D)
                matrices[0][verCell][horCell].value = value_P

                if (value_P == up_P):
                    matrices[0][verCell][horCell].source.append(
                        "0," + str(verCell - 1) + "," + str(horCell) + ",UP")
                if (value_P == up_D):
                    matrices[0][verCell][horCell].source.append(
                        "2," + str(verCell - 1) + "," + str(horCell) + ",UP")

                valueDiagD = matrices[2][verCell - 1][horCell - 1].value + int(self.scoringMatrix[self.sequenceA[horCell - 1]][self.sequenceB[verCell - 1]])
                valueD = max(valueDiagD, value_P, value_Q)
                matrices[2][verCell][horCell].value = valueD

                if(valueD == value_P):
                    matrices[2][verCell][horCell].source.append("0," + str(verCell) + "," + str(horCell) + ",P")
                if(valueD == value_Q):
                    matrices[2][verCell][horCell].source.append("1," + str(verCell) + "," + str(horCell) + ",Q")
                if(valueD == valueDiagD):
                    matrices[2][verCell][horCell].source.append("2," + str(verCell - 1) + "," + str(horCell - 1) + ",DIAGONAL")

                horCell+= 1
            verCell += 1
            horCell = 1

        self.finalAlignScore = matrices[2][len(self.sequenceB)][len(self.sequenceA)].value

        return matrices
    """
    Traceback for Gotoh algorithm. If value comes from upper cell then add '-', if left then '-'
    and match/mismatch for diagonal
    """
    def Traceback(self):
        matrices = self.filledMatrices
        seqA = self.sequenceA
        seqB = self.sequenceB
        
        column = len(seqA)
        Row = len(seqB)
        SequenceA = len(seqA) - 1
        SequenceB = len(seqB) - 1
        constSeqA = []
        constSeqB = []
        loc = 2

        while Row != 0 or column !=0 or loc != 2:
            source = matrices[loc][Row][column].source
            comp = source[random.randint(0, len(source) - 1)].split(',')
            loc = int(comp[0])
            Row = int(comp[1])
            column = int(comp[2])

            if comp[3] == "UP":
                constSeqA.insert(0, '-')
                constSeqB.insert(0, seqB[SequenceB])
                SequenceB -= 1
            elif comp[3] == "LEFT":
                constSeqA.insert(0, seqA[SequenceA])
                constSeqB.insert(0, '-')
                SequenceA -= 1
            elif comp[3] == "DIAGONAL":
                constSeqA.insert(0, seqA[SequenceA])
                constSeqB.insert(0, seqB[SequenceB])
                SequenceA -= 1
                SequenceB -= 1

        #join them together
        alignedSeqA = ''.join(constSeqA)
        alignedSeqB = ''.join(constSeqB)

        return [alignedSeqA, alignedSeqB]
    
    # for Error output
    def Errors(self, error):
        if (isinstance(error, list)):
            for err in error:
                print(err)
        elif (isinstance(error, str)):
            print(error)

    # # for pairwise alignment output
    def NWAlign(self, Result):
        if (not isinstance(Result[0], list)):
            resultA = self.output(Result[0], 80)
            resultB = self.output(Result[1], 80)
            seperator = self.Seperator(Result[0], Result[1])
            resultFormat = self.output(seperator, 80)

            for i in range(0, len(resultA)):
                a = i + 1
                print("Aligned Sub-Part No:", a, "\n")
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


def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastaFile", type=str)
    parser.add_argument("-gap", "--gapPenality", type=int, default=-4)
    parser.add_argument("-ge", "--gapExtendCost", type=int, default=-1)
    parser.add_argument("-matrix", "--scoringMatrix", type=str, default="pam250")
    args = parser.parse_args()

    readSequences = FileReader.fastaFile(args.fastaFile)
    seqA = readSequences[0]
    seqB = readSequences[1]

    align = Gotoh()
    align.sequenceA = seqA
    align.sequenceB = seqB
    align.gapPenality = args.gapPenality
    align.IP_Matrix = args.scoringMatrix
    align.gapExtended = args.gapExtendCost
    errors = align.errorCheck()
    if (errors == None):
        Result = align.gotoh()
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("\t" + "\t" + "\t" + "Input Sequences for Gotoh: " + "\n")
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("Sequence A: " + "\n" + "\n" + seqA + "\n")
        print("Sequence B: " + "\n" + "\n" + seqB + "\n")
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("Gap Penality: " + str(args.gapPenality) + "\n")
        print("Scoring Matrix : " + align.IP_Matrix + "\n")
        print("Extended Gap : " + str(align.gapExtended) + "\n")
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("\t" + "\t" + "\t" + "Output sequence is divided into 80 letters" + "\n")
        align.NWAlign(Result)
        print("_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ")
        print("Alignment Score: " + str(align.finalAlignScore))
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    else:
        print("\n" + "_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  Error  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _" + "\n")
        align.Errors(errors)


if __name__ == "__main__":
    Main()
