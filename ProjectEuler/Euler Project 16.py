num_str = map (int, str (2**1000))

from operator import add

sum = reduce (add, num_str)

print sum
