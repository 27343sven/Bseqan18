import numpy as np

class Aligner:

    def __init__(self, seq1="", seq2=""):
        self.seq1 = seq1
        self.seq2 = seq2
        self.penalties = score_kimura
        self.scoreMatrix = np.zeros((len(self.seq1), len(self.seq2)), dtype=np.int)
        self.directionMatrix = np.zeros((len(self.seq1), len(self.seq2)), dtype=np.int)

    def needleman_wunch(self, seq1="", seq2=""):
        self.seq1 = seq1 if seq1 else self.seq1
        self.seq2 = seq2 if seq2 else self.seq2
        self.__checkSequence()
        self.seq1 = "*" + self.seq1
        self.seq2 = "*" + self.seq2
        self.__fill_matrices()
        self.__calulate_scores()
        answer = self.__traceback()
        self.seq1 = self.seq1[1:]
        self.seq2 = self.seq2[1:]
        return answer

    def smith_waterman(self, seq1="", seq2=""):
        self.seq1 = seq1 if seq1 else self.seq1
        self.seq2 = seq2 if seq2 else self.seq2
        self.__checkSequence()
        self.seq1 = "*" + self.seq1
        self.seq2 = "*" + self.seq2
        self.__fill_matrices(waterman=True)
        self.__calulate_scores(waterman=True)
        answer = self.__traceback(waterman=True)
        self.seq1 = self.seq1[1:]
        self.seq2 = self.seq2[1:]
        return answer

    def __calulate_scores(self, waterman=False):
        for i in range(1, self.scoreMatrix.shape[0]):
            for j in range(1, self.scoreMatrix.shape[1]):
                x = self.scoreMatrix.shape[1]
                values = np.take(self.scoreMatrix, [i*x + j-1, (i-1)*x + j-1, (i-1)*x + j])
                modifiers = np.array([-5, self.penalties[self.seq1[i] + self.seq2[j]], -5])
                values += modifiers
                if waterman and values.max() < 0:
                    self.directionMatrix[i, j] = 0
                    self.scoreMatrix[i, j] = 0
                else:
                    self.directionMatrix[i, j] = int("".join(list(np.where(values == np.max(values), "1", "0"))), 2)
                    self.scoreMatrix[i, j] = np.max(values)

    def __traceback(self, waterman=False):
        x = self.scoreMatrix.shape[1]
        _seq1, _seq2, i, j = "", "", self.scoreMatrix.shape[0]-1, self.scoreMatrix.shape[1]-1
        if waterman:
            maximum = self.scoreMatrix.argmax()
            i, j = maximum // x, maximum % x
        while self.scoreMatrix[i, j] != 0 if waterman else (i, j) != (0, 0):
            directions = np.array(list(bin(self.directionMatrix[i, j])[2:].zfill(3)), dtype=np.int)
            values = np.take(self.scoreMatrix, [max(0, i * x + j - 1),
                                                max(0, (i - 1) * x + j - 1),
                                                max(0, (i - 1) * x + j)])
            values = np.ma.masked_array(values, np.logical_not(directions))
            np.ma.set_fill_value(values, min(values) - 1)
            new_dir = np.ma.argmax(values)
            _seq1 += "_" if new_dir == 0 else self.seq1[i]
            _seq2 += "_" if new_dir == 2 else self.seq2[j]
            i = i - 0 if new_dir == 0 else i - 1
            j = j - 0 if new_dir == 2 else j - 1
        answer = AlignementResult(self.seq1[1:], self.seq2[1:], _seq1[::-1], _seq2[::-1], self.scoreMatrix, self.directionMatrix)
        return answer

    def __fill_matrices(self, waterman=False):
        self.scoreMatrix = np.zeros((len(self.seq1), len(self.seq2)), dtype=np.int)
        self.directionMatrix = np.zeros((len(self.seq1), len(self.seq2)), dtype=np.int)
        if waterman:
            self.scoreMatrix[0] = np.full(self.scoreMatrix.shape[1], 0)
            self.scoreMatrix[:, 0] = np.full(self.scoreMatrix.shape[0], 0)
        else:
            self.scoreMatrix[0] = np.arange(0, self.scoreMatrix.shape[1] * self.penalties["*"], self.penalties["*"])
            self.scoreMatrix[:, 0] = np.arange(0, self.scoreMatrix.shape[0] * self.penalties["*"], self.penalties["*"])
        self.directionMatrix[0, 1:] = np.full(self.directionMatrix.shape[1] - 1, int('100', 2))
        self.directionMatrix[1:, 0] = np.full(self.directionMatrix.shape[0] - 1, int('001', 2))

    def __checkSequence(self):
        self.seq1 = self.seq1.upper()
        self.seq2 = self.seq2.upper()

        if self.seq1 == "" or self.seq2 == "":
            raise SequenceException("One or both sequences are empty.")
        if self.seq1.strip("ATCG") != "" or self.seq2.strip("ATCG") != "":
            raise SequenceException("One or both of the sequences contain illegal characters.")

class AlignementResult:

    def __init__(self, seq1, seq2, al_seq1, al_seq2, scoreMatrix, directionMatrix):
        self.seq1 = seq1
        self.seq2 = seq2
        self.al_seq1 = al_seq1
        self.al_seq2 = al_seq2
        self.scoreMatrix = scoreMatrix
        self.directionmatrix = directionMatrix
        self.al_score = {'match': 5, 'mis': -4, 'gap': -10, 'ext': -1}

    def show_scorematrix(self):
        print(self.scoreMatrix)

    def show_directionmatrix(self):
        print(self.directionmatrix)

    def show_alignment(self):
        print(self.al_seq1)
        print(self.__make_alignment_symbols())
        print(self.al_seq2)
        print("score: {}".format(self.__calc_alignmentscore()))

    def save_alignment(self, filename):
            with open(filename if filename[-4:] == ".txt" else filename + ".txt", "w") as file:
                file.write(self.al_seq1 + "\n")
                file.write(self.__make_alignment_symbols() + "\n")
                file.write(self.al_seq2 + "\n")
                file.write("score: {}\n".format(self.__calc_alignmentscore()))

    def save_traceback(self, filename):
        self.seq1 = "*" + self.seq1
        with open(filename if filename[-4:] == ".csv" else filename + ".csv", "w") as file:
            file.write(";*;{}\n".format(";".join(self.seq2)))
            for i in range(len(self.seq1)):
                dir = ["[{1}] {0}".format(x, self.__make_direction(self.directionmatrix[i, j])) for j, x in enumerate(list(self.scoreMatrix[i]))]
                file.write("{};{}\n".format(self.seq1[i], ";".join(dir)))
        self.seq1 = self.seq1[1:]

    def __make_direction(self, i):
        direction = ""
        binary = list(bin(i)[2:].zfill(3))
        if binary[0] == "1":
            direction += "-"
        if binary[1] == "1":
            direction += "\\"
        if binary[2] == "1":
            direction += "|"
        if binary == ["0", "0", "0"]:
            direction += "*"
        return direction

    def __make_alignment_symbols(self):
        symbols = ""
        for i in range(len(self.al_seq1)):
            if self.al_seq1[i] == self.al_seq2[i]:
                symbols += "|"
            elif self.al_seq1[i] == "_" or self.al_seq2[i] == "_":
                symbols += " "
            else:
                symbols += "*"
        return symbols

    def __calc_alignmentscore(self):
        score, gap = 0, False
        for i in range(len(self.al_seq1)):
            if self.al_seq1[i] == self.al_seq2[i]:
                score += self.al_score['match']
                if gap:
                    gap = False
            elif self.al_seq1[i] == "_" or self.al_seq2 == "_":
                if gap:
                    score += self.al_score['ext']
                else:
                    gap = True
                    score += self.al_score['gap']
            else:
                score += self.al_score['mis']
                if gap:
                    gap = False
        return score





class SequenceException(Exception):

    def __init__(self, message):
        super().__init__(message)


score_jukes_cantor = {"AA": 5, "AT": -4, "AC": -4, "AG": -4,
                      "TT": 5, "TA": -4, "TC": -4, "TG": -4,
                      "CC": 5, "CA": -4, "CT": -4, "CG": -4,
                      "GG": 5, "GA": -4, "GT": -4, "GC": -4,
                      "A*": -5, "T*": -5, "C*": -5, "G*": -5,
                      "*A": -5, "*T": -5, "*C": -5, "*G": -5,
                      "**": 0, "*": -5}

score_kimura = {"AA": 5, "AT": -4, "AC": -4, "AG": 0,
                "TT": 5, "TA": -4, "TC": 0, "TG": -4,
                "CC": 5, "CA": -4, "CT": 0, "CG": -4,
                "GG": 5, "GA": 0, "GT": -4, "GC": -4,
                "A*": -5, "T*": -5, "C*": -5, "G*": -5,
                "*A": -5, "*T": -5, "*C": -5, "*G": -5,
                "**": 0, "*": -5}
