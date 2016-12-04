import re
import ast
import sys

with open(sys.argv[1]) as f:
    total = 0
    count = 0
    for line in f:
        total += float(line)
        count += 1

average = (total/count)/1000000
print(average)
print(count)