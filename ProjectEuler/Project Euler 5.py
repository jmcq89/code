max = 20*19*18*17*16*15*14*13*12*11

def isMult(n):
    if n%19==0 and n%18==0 and n%17==0 and n%16==0 and n%15==0 and \
    n%14==0 and n%13==0 and n%12==0 and n%11==0:
        return True

def lowestMult():
    n=20
    while n<max:
        if isMult(n):
            print n
            break
        n+=20

lowestMult()

print isMult(232792560)
