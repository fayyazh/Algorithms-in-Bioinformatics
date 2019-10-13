import FileReader

class Sum_Of_Pairs:
    def __init__(self, alignment, Scoring_Function):
        self.alignment = alignment
        self.Scoring_Function = Scoring_Function
        self.errors = []
        self.file = []

    def ComputeScore(self, alignment = None):
        errors = []

        try:
            self.file = FileReader.ScoringMatrix(self, self.Scoring_Function)
        except FileNotFoundError:
            errors.append("file is not found: " + self.Scoring_Function)
        except IOError:
            errors.append("file cannot be opened: " + self.Scoring_Function)

        if (len(errors) == 0):
            self.errors = None
            self.file = FileReader.ScoringMatrix(self, self.Scoring_Function)
        else:
            self.errors = "\n".join(errors)

        SOP = 0

        i = 1
        for seqA in alignment:
            for seqB in alignment[i:]:
                Score = 0

                for j in range(len(seqA)):
                    SequenceA = seqA[j]
                    SequenceB = seqB[j]

                    SequenceA = SequenceA.replace('-', '*')
                    SequenceB = SequenceB.replace('-', '*')

                    Score += int(self.file[SequenceA] [SequenceB])

                SOP += Score

            i += 1

        return SOP
