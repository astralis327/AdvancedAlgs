#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
    # handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, sequences, table, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        results = []

        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):

                align_lengthi, align_lengthj = align_length, align_length
                horiz_string = "_" + sequences[i]       #insert a blank character for the base cases
                vert_string = "_" + sequences[j]
                if (align_lengthi > len(horiz_string)): align_lengthi = len(horiz_string)   #if the align length is bigger than the string length

                if align_lengthj > len(vert_string):
                    align_lengthj = len(vert_string)
                if (j < i):
                    s = {}

                else:
                    if not banded:  #normal unrestricted alg
                        score, alignment1, alignment2 = \
                            self.unrestricted(horiz_string, vert_string, align_lengthi, align_lengthj, False)

                    else:

                        # if doing banded alignment, diagonal will end long before the string does if matching the 1000 char string
                        # against "polynomial" or "exponential."
                        if banded:
                            if align_lengthj > align_lengthi + 3:
                                align_lengthj = align_lengthi

                        score, alignment1, alignment2 = \
                            self.banded_alignment(horiz_string, vert_string, align_lengthj)
                    ###################################################################################################
                    # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
                    # 					score = i+j;
                    # 					alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                    # 						len(sequences[i]), align_length, ',BANDED' if banded else '')
                    # 					alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                    # 						len(sequences[j]), align_length, ',BANDED' if banded else '')
                    ###################################################################################################
                    s = {'align_cost': score, 'seqi_first100': alignment1[:100], 'seqj_first100': alignment2[:100]}
                    table.item(i, j).setText('{}'.format(int(score) if score != math.inf else score))
                    table.update()
                jresults.append(s)
            results.append(jresults)
        return results

    # my understanding is that gaps are taking scores from the left or from above, and substitutions
    # /single char mismatches take from the diagonal.

    # Accidentally made a column by column matrix instead of a row by row matrix, so basically a transposed\
    # version of what it's supposed to be. Answers should still be the same, however.

    def banded_alignment(self, string_horiz, string_vert, len_j):
        num_cols = 2 * 3 + 1       #number of columns in the transformed matrix
        score = []  # score list of rows
        pointers = []  # pointers matrix of the directions the next best alignment is in
        # num rows one more integer than the length of the string to account for the base cases
        d = 3
        for row_num in range(len_j):
            counter = 0                 #counter to tell whether the null corners are passed
            row = []
            pointer_row = []
            transform = -1 * row_num + d    #cell transformation

            score.append(row)
            pointers.append(pointer_row)
            for col in range(num_cols):
                colNorm = col - transform   #normal column number in untransformed matrix
                if counter < transform and row_num >= 0 \
                        or counter > -1 * transform and row_num >= len_j - 1 - 2:   #if the corners of the matrix have been reached
                    score[row_num].insert(col, math.inf)  # empty cells in the corners
                    pointers[row_num].insert(col, "None")
                elif row_num == 0 and counter >= transform:
                    score[row_num].insert(col, col * 5 - d * 5)     #base case at the top row
                    pointers[row_num].insert(col, "insert")

                elif colNorm < len(string_horiz):                   #check to make sure the col num isn't out of bounds of its string
                    if string_horiz[colNorm] == string_vert[row_num]:
                        cost = score[row_num - 1][col] - 3  # match, diagonal cel in normal matrix
                        pointers[row_num].insert(col, 'match')
                    else:
                        pointers, cost = self.find_minimum(score, col, row_num, pointers, True)     #compute cost if not a match
                    score[row_num].insert(col, cost)

                counter += 1
        horiz_align, vert_align = self.create_alignment(pointers, string_horiz, string_vert, True) #get the string alignment
        return score[row_num][len(score[row_num])-1], horiz_align, vert_align

    def unrestricted(self, string_horiz, string_vert, len_i, len_j, banded):

        score = []  # score list of columns
        pointers = []  # pointers matrix of the directions the next best alignment is in
        rows, cols = len_j, len_i  # num rows one more integer than the length of the string to account for the base cases

        for k in range(rows):  # filling out the matrix column by column
            for l in range(cols):
                row = []
                pointer_row = []
                cost = 0
                if k == 0:  # first row
                    row.append(l * 5)  # base case, col num multiplied by 5 to mimic insertions
                    pointer_row.append("insert")
                elif l == 0:  # first col
                    score[0].insert(k, k * 5)      #fill first col up with base cases
                    pointers[0].insert(k, "del")
                else:
                    if string_horiz[l] == string_vert[k]:
                        cost = score[l - 1][k - 1] - 3  # match
                        pointers[l].insert(k, "match")
                    else:
                        pointers, cost = self.find_minimum(score, l, k, pointers, False)    #calculate the cost if not a match
                    score[l].insert(k, cost)

                if len(row) > 0:
                    score.append(row)
                    pointers.append(pointer_row)
        horiz_align, vert_align = self.create_alignment(pointers, string_horiz, string_vert, False)
        return score[cols - 1][rows - 1], horiz_align, vert_align

    def find_minimum(self, score, col, row, pointers, banded):

        if banded:
            if col == 6: choice1 = math.inf		#edge case, right side of matrix
            else:
                choice1 = score[row - 1][col + 1] + 5  # top cell, deletion (adding gap in horiz string)

            if col == 0: choice2 = math.inf						#edge case, left side of matrix
            else:
                choice2 = score[row][col - 1] + 5  # left cell, insertion (adding gap in vert string
            choice3 = score[row-1 ][col] + 1  # diagonal, subsitution
            dim1, dim2 = row, col
        else:
            dim1, dim2 = col, row
            choice2 = score[col - 1][row] + 5  # left cell, insertion (adding gap in vert string)
            choice1 = score[col][row - 1] + 5  # top cell, deletion (adding gap in horiz string)
            choice3 = score[col - 1][row - 1] + 1  # diagonal, subsitution

        minimum = choice2  # set the left cell as the default choice
        pointers[dim1].insert(dim2, "insert")
        if choice1 < minimum:
            minimum = choice1
            pointers[dim1][dim2] = "del"
        if choice3 < minimum:
            minimum = choice3
            pointers[dim1][dim2] = "sub"
        return pointers, minimum

    def create_alignment(self, pointers, horiz_string, vert_string, banded):
        if banded:
            row = len(pointers) - 1             #rows and columns for the banded and unrestricted algs are different because
            col = len(pointers[row]) -1         #the unrestricted matrix was set up as an array of columns and the banded matrix
        else:                                   #is a list of rows, but both the banded and unrestricted algs return the same answer for
            col = len(pointers) - 1             #basic test cases and I dont have time to change it
            row = len(pointers[col]) - 1


        while row > 0 and col > 0:

            if banded:
                colNorm = col - (-1 * row + 3)  #untransformed col num
                element = pointers[row][col]    #get the direction to go in (left for insert, up and right for del, and up for match or sub)
                # insert - in horiz string
                if element == "insert":
                    col -= 1
                if element == "del":
                    row -= 1
                    col += 1
                elif element == "sub" or element == "match":
                    row -= 1

            else:
                colNorm = col
                element = pointers[col][row]
             # insert - in horiz string
                if element == "insert": #left for insert, up for del, and up and right for sub or match
                    col -= 1
                if element == "del":
                    row -= 1
                elif element == "sub" or element == "match":
                    col -= 1
                    row -= 1
            if element == "del":
                horiz_string = horiz_string[:(colNorm+1)] + "-" + horiz_string[(colNorm+1):] #insert a gap in the horizontal string if there's a deletion
            elif element == "insert":
                vert_string = vert_string[:(row+1)] + "-" + vert_string[(row+1):]   #insert a gap in vert string if there's an insertion

        horiz_string = horiz_string[1:]     #delete the underscores put in at the beginning
        vert_string = vert_string[1:]
        return horiz_string, vert_string
