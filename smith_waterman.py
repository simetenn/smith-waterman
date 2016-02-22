from Bio import SeqIO, SubsMat
import numpy as np
from Bio.SubsMat.MatrixInfo import blosum62
import sys






def loadScoringMatrix(scoring_matrix_file):
    """
    Load the scoring matrix from file and into a double dictionary

    The format must be as the example blossum62.txt file
    """
    tmp_scoring_matrix = np.loadtxt(scoring_matrix_file, str)

    symbols = tmp_scoring_matrix[0, 1:]
    tmp_scoring_matrix = tmp_scoring_matrix[1:, 1:]

    scoring_matrix = {}
    i = 0
    for symbol_i in symbols:
        scoring_matrix[symbol_i] = {}
        j = 0
        for symbol_j in symbols:
            scoring_matrix[symbol_i][symbol_j] = float(tmp_scoring_matrix[i][j])
            j += 1

        i += 1

    return scoring_matrix


def test_similarity(a, b):
    if a == b:
        return 2
    else:
        return -1

def smithWatherman(sequence1, sequence2, scoring_matrix_file, W, k=None):
    scoring_matrix = loadScoringMatrix(scoring_matrix_file)

    m = len(sequence1)
    n = len(sequence2)

    H = np.zeros((m + 1, n + 1))

    for i in xrange(1, m + 1):
        for j in xrange(1, n + 1):
            # print "-------------------------"
            # diagonal_score = H[i-1, j-1] + scoring_matrix[sequence1[i - 1]][sequence2[j - 1]]
            diagonal_score = H[i-1, j-1] + test_similarity(sequence1[i - 1], sequence2[j - 1])

            if k is None:
                k = np.arange(1, i + 1)
            left_score = (H[i - k, j] + W(k)).max()

            if k is None:
                k = np.arange(1, j + 1)
            up_score = (H[i, j - k] + W(k)).max()

            # print sequence1[i - 1], sequence2[j - 1]
            # print diagonal_score, left_score, up_score

            H[i, j] = max(0, diagonal_score, up_score, left_score)
            # print H[i, j]
        # sys.exit(0)
    return H.T


def traceback(sequence1, sequence2, H):
    x, y = np.unravel_index(np.argmax(H), H.shape)

    aligned_sequence1 = ""
    aligned_sequence2 = ""

    print x, y, H[x, y]


def trace(x, y, H, sequence1, sequence2, aligned_sequence1, aligned_sequence2):

        if H[x - 1, y - 1] >= max(H[x - 1, x], H[x, y - 1]):
            x -= 1
            y -= 1

        if H[x - 1, x]:
            x -= 1

        if H[x, y - 1]:
            y -= 1

if __name__ == "__main__":


    filename1 = "NP_001264426.fasta"
    filename2 = "Q6NS38.1.fasta"

    scoring_matrix_file = "blossum62.txt"

    sequence1 = SeqIO.read(open(filename1), 'fasta').seq
    sequence2 = SeqIO.read(open(filename2), 'fasta').seq

    sequence1 = "ACACACTA"
    sequence2 = "AGCACACA"

    sequence1 = "CTCATGC"
    sequence2 = "ACAATCG"


    def W(length):
        return -(11 + length)

    def test_W(length):
        return -1

    H = smithWatherman(sequence1, sequence2, scoring_matrix_file, test_W, k=1)
    traceback(sequence1, sequence2, H)
