file = open('Euler13number.txt', 'r')
print str(sum([int(line) for line in file]))[:10]
file.close()
