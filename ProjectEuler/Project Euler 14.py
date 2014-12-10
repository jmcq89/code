def lengthCollatz(n):
    i=1
    while n>1:
        i+=1
        if n%2==1:
            n=3*n+1
        else:
            n=n/2
    return i

def longestCollatz(max):
    length=0
    num=0
    for i in range(max-1,500000,-2):
        if lengthCollatz(i)>length:
            length=lengthCollatz(i)
            num=i
    print length
    print num

longestCollatz(1000000)
