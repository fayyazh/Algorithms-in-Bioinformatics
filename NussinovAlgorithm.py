import FileReader
import argparse
import math

class MatrixConst:
    def __init__(self, val=0):
        self.value = val
        self.source = []

    def __repr__(self):
        return str(self.source)

class NussinovAlgorithm:

    def __init__(self, sequence, Scoring_Function):
        self.sequence = sequence.upper()
        self.location = self.Location()
        self.Scoring_Function = Scoring_Function
        self.scoringMatrix = []
        self.endMatrix = None
        self.finalAlignment = None


    def errorCheck(self):
        errors = []

        nt = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
        for seq in self.sequence:
            for letter in seq:
                if (letter not in nt):
                    errors.append("Sequence contains invalid letter: " + letter)

        try:
            self.scoringMatrix = FileReader.ScoringMatrix(self, self.Scoring_Function)
        except FileNotFoundError:
            errors.append("file is not found: " + self.Scoring_Function)
        except IOError:
            errors.append("file is not readable: " + self.Scoring_Function)

        if (len(errors) == 0):
            self.errors = None
        else:
            self.errors = "\n".join(errors)

        return errors


    def Location(self, sequence = None):
        if sequence is None: sequence = self.sequence

        location = {}
        pos = 0
        for char in sequence:
            if char not in location:
                location[char] = [pos]
            else:
                location[char].append(pos)
            pos += 1

        return location
    
    def Nussinov(self):
        self.endMatrix = self.Matrix_Filling()
        self.finalAlignment = self.Traceback()
        return self.finalAlignment
    
    def Matrix_Initialization(self, sequence = None):
        if sequence is None: sequence = self.sequence

        # constructing matrix
        matrix = [[MatrixConst(None) for i in range(len(sequence))] for j in range(len(sequence))]
        
        matrix[0][0].value = 0
        i = 1
        j = 1
        while i <= len(sequence) - 1:
            matrix[i][j].value = 0
            matrix[i][j - 1].value = 0
            i +=1
            j +=1

        return matrix
    

    def Matrix_Filling(self, sequence = None):
        if sequence is None: sequence = self.sequence

        matrix = self.Matrix_Initialization()

        Current_Row = 0
        Current_Col = 1
        Start_Col = 1

        while Start_Col < len(sequence):
            #fill diagonal
            while Current_Col < len(sequence):
                Left = matrix[Current_Row][Current_Col - 1].value
                Down = matrix[Current_Row + 1][Current_Col].value
                Diagonal = matrix[Current_Row + 1][Current_Col - 1].value + int(self.scoringMatrix[sequence[Current_Row]][sequence[Current_Col]])


                complementary_Loc = self.Complementary_Loc(sequence[Current_Col])
                Value = float("-inf")
                Source = []
                for pos in complementary_Loc:
                    if(pos > Current_Row and pos < Current_Col):
                        val = matrix[Current_Row][pos].value + matrix[pos + 1][Current_Col].value
                        if(val > Value):
                            Value = int(val)
                            Source = [str(Current_Row) + "," + str(pos), str(pos + 1) + "," + str(Current_Col)]
                        elif(val == Value):
                            Source.append(str(Current_Row) + "," + str(pos))
                            Source.append(str(pos + 1) + "," + str(Current_Col))
                seperation = Value

                #determine max
                value = max(Left, Down, Diagonal, seperation)
                matrix[Current_Row][Current_Col].value = value

                if(value == Diagonal):
                    matrix[Current_Row][Current_Col].source.append(str(Current_Row + 1) + "," + str(Current_Col - 1) + ",DIAGONAL")
                if(value == Left):
                    matrix[Current_Row][Current_Col].source.append(str(Current_Row) + "," + str(Current_Col - 1))
                if(value == Down):
                    matrix[Current_Row][Current_Col].source.append(str(Current_Row + 1) + "," + str(Current_Col))
                if(value == seperation):
                    for src in Source:
                        matrix[Current_Row][Current_Col].source.append(src)

                Current_Row += 1
                Current_Col += 1
            Start_Col += 1
            Current_Row = 0
            Current_Col = int(Start_Col)

        return matrix

    def Complementary_Loc(self, base):
        try:
            return self.location[self.Complementary_Base(base)]
        except KeyError:
            return []
        
    def Complementary_Base(self, base):
        if(base == "C"):
            return "G"
        elif(base == "G"):
            return "C"
        elif(base == "A"):
            return "U"
        elif(base == "U"):
            return "A"


   # Unable to implement the traceback correctly
    def Traceback(self, matrix = None):
        if matrix is None: matrix = self.endMatrix
        
        Output = []
        Current_Row = 0
        Current_Col = len(self.sequence) - 1
        
        while(True):
            if(len(matrix[Current_Row][Current_Col].source) > 0):
                comp = matrix[Current_Row][Current_Col].source[0].split(",")
                if(len(comp) == 3 and comp[2] == "DIAGONAL"):
                    Output.append([Current_Row, Current_Col])
                Current_Row = int(comp[0])
                Current_Col = int(comp[1])
            else:
                break

        return Output

    def PrintNussinov(self, sequence, Output):
        comma = ',' * len(sequence)
        commaList = list(comma)
        for pair in Output:
            commaList[pair[0]] = '('
            commaList[pair[1]] = ')'
        comma = ''.join(commaList)

        Sequence = self.output(sequence, 80)
        Seperator = self.output(comma, 80)

        for i in range(len(Sequence)):
            print(Sequence[i])
            print(Seperator[i])

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
    parser.add_argument("-s", "--Scoring_Function", type=str, default="pam250")
    Pars = parser.parse_args()

    error = []
    sequences = FileReader.fastaFile(Pars.fastaFile)

    Nuss = NussinovAlgorithm(sequences[0], Pars.Scoring_Function)
    error = Nuss.errorCheck()
    if len(error) == 0:
        Output = Nuss.Nussinov()

        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("NEEDLEMANWUNSCH ALGORITHM".center(80, ' ') + "\n")
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("INPUT SEQUENCES".center(80, ' ') + "\n")
        print("Sequence: " + Nuss.sequence)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("SCORING FUNCTION".center(80, ' '))
        print("Scoring Matrix: " + Nuss.Scoring_Function)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print("OUTPUT".center(80, ' ') + "\n")
        Nuss.PrintNussinov(Nuss.sequence, Output)
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    else:
        print("Errors found:")
        Nuss.Errors(error)


if __name__ == "__main__":
    Main()













