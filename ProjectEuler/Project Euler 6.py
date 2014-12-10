def isPrime(n):
    n*=1.0
    if n%2==0 and n!=2 or n%3==0 and n!=3:
        return False
    for b in range(1,int((n**0.5+1)/6.0+1)):
        if n%(6*b-1)==0:
            return False
        if n %(6*b+1)==0:
           return False
    return True

def nthPrime(n):
    nth=1 # we know 2 is prime
    num=3 # 
    while nth<n+1:
        if isPrime(num):
            nth+=1
            if nth==n:
                print num
        num+=2

nthPrime(10001)
