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

def sumPrimesUpTo(max):
    count = 1
    sum = 2
    num = 3
    while num < max+1:
        if isPrime(num):
            sum+=num
            count+=1
        num+=2
    print sum
    
sumPrimesUpTo(2000000)

