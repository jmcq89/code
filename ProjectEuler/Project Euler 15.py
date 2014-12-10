def nthPascal(n):
    print 1
    print n
    num=(n*(n-1))/2
    print num
    i=2
    while num>1:
        num=(num*(n-i))/(i+1)
        print num
        i+=1

nthPascal(40)
        
