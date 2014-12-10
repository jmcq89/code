## Project Euler Problem 21
## Find the sum of Amicable pairs under 10000

from math import sqrt

def factors(n):
    fact=[1]
    check=2
    rootn=sqrt(n)
    while check<rootn:
    	if n%check==0:
    		fact.append(check)
    		fact.append(n/check)
    	check+=1
    if rootn==check:
    	fact.append(check)
    	fact.sort()
    return sorted(fact)

def findAmicablePairs(n):
    i=200
    pairs=[]
    while i<n+1:
        a=sum(factors(i)) ## 284
        b=sum(factors(a)) ## 220
        if b==i and a!=i:
            pairs.append(i)
        i+=1
    return pairs

print sum(findAmicablePairs(10000))
