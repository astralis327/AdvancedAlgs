

import sys

sequence = [3,6,8,-6,30,-30,48]
index = 0
sum = 0
maximum= []
for i in range(len(sequence)):
    maximum.append(sum+sequence[i])
    if maximum[i] > sum:
        index = i
        sum = sequence[i] + sum
    elif maximum[i] <= 0:
        maximum[i] = 0

print(index)