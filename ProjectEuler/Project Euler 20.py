## Project Euler Problem 20:
## http://projecteuler.net/index.php?section=problems&id=20
## Find the sum of the digits of 100!

import math
from operator import add

num=math.factorial(100)
num_str=map(int,str(num))
sum=reduce(add,num_str)
print sum 
