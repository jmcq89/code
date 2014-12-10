from math import sqrt

def factors(n):
    fact=[1,n]
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
    return fact

def triangleFactors(max):
    n=6
    triangle=21
    while len(factors(triangle))<max+1:
        triangle=(n*(n+1))/2
        n+=1
        if len(factors(triangle))>max:
            print triangle
            

triangleFactors(500)

