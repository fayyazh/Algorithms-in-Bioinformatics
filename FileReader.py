def ScoringMatrix(self, fileName):
        file = open(fileName + ".txt", "r")
        linesRead = file.readlines()
        file.close()

        Dictionary = {}
        sequence = []
        for lineNo in linesRead:
            if lineNo[0] != '#':
                if lineNo[0] == ' ':
                    string = lineNo.strip(' ').replace('\n', '')
                    sequence = string.split('  ')
                else:
                    string = lineNo.strip(' ').replace('\n', '').replace('  ', ' ')
                    row = string.split(' ')
                    rDict = {}
                    i = 0
                    for value in row[1:]:
                        rDict[sequence[i]] = int(value)
                        i += 1
                    Dictionary[row[0]] = rDict
        return Dictionary

def fastaFile(filename):
        string = ""
        sequences = []
        file = open(filename + ".fasta", "r")
        linesRead = file.readlines()
        file.close()
        for line in linesRead:
                if(line[0] == ">"):
                        sequences.append(string)
                        string = ""
                else:
                        string += line.replace("\n", "")
        sequences.append(string)
        del sequences[0]
        return sequences

