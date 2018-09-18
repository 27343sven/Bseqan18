import numpy as np

class Aligner:

    def __init__(self, seq1="", seq2=""):
        self.seq1 = seq1
        self.seq2 = seq2
        self.penalties = score_kimura
        self.scoreMatrix = np.zeros((len(self.seq1), len(self.seq2)), dtype=np.int)
        self.directionMatrix = np.zeros((len(self.seq1), len(self.seq2)), dtype=np.int)

    def smith_waterman(self, seq1="", seq2=""):
        self.seq1 = seq1 if seq1 else self.seq1
        self.seq2 = seq2 if seq2 else self.seq2
        self.__checkSequence()
        self.seq1 = "*" + self.seq1
        self.seq2 = "*" + self.seq2
        self.__fill_matrices()
        self.__calulate_scores()
        print(self.scoreMatrix)
        print(self.directionMatrix)
        self.__traceback()

    def __calulate_scores(self):
        for i in range(1, self.scoreMatrix.shape[0]):
            for j in range(1, self.scoreMatrix.shape[1]):
                x = self.scoreMatrix.shape[1]
                values = np.take(self.scoreMatrix, [i*x + j-1, (i-1)*x + j-1, (i-1)*x + j])
                modifiers = np.array([-5, self.penalties[self.seq1[i] + self.seq2[j]], -5])
                values += modifiers
                self.directionMatrix[i, j] = int("".join(list(np.where(values == np.max(values), "1", "0"))), 2)
                self.scoreMatrix[i, j] = np.max(values)

    def __traceback(self):
        x = self.scoreMatrix.shape[1]
        _seq1, _seq2, i, j = "", "", self.scoreMatrix.shape[0]-1, self.scoreMatrix.shape[1]-1
        while (i, j) != (0, 0):
            print("===({}, {})===".format(i, j))
            directions = np.array(list(bin(self.directionMatrix[i, j])[2:].zfill(3)), dtype=np.int)
            values = np.take(self.scoreMatrix, [max(0, i * x + j - 1),
                                                max(0, (i - 1) * x + j - 1),
                                                max(0, (i - 1) * x + j)])
            print(directions)
            print(values)
            values = np.ma.masked_array(values, np.logical_not(directions[directions]))
            np.ma.set_fill_value(values, min(values) - 1)
            new_dir = np.ma.argmax(values)
            print(new_dir)
            _seq1 += "_" if new_dir == 0 else self.seq1[i]
            _seq2 += "_" if new_dir == 2 else self.seq2[j]
            i = i - 0 if new_dir == 0 else i - 1
            j = j - 0 if new_dir == 2 else j - 1
            print(_seq1[::-1])
            print(_seq2[::-1])





    def __fill_matrices(self):
        self.scoreMatrix = np.zeros((len(self.seq1), len(self.seq2)), dtype=np.int)
        self.directionMatrix = np.zeros((len(self.seq1), len(self.seq2)), dtype=np.int)
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
