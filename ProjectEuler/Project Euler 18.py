def _reduce_triangle(to_reduce):
    last_row=to_reduce[-1]
    for index in xrange(len(to_reduce)-1):
        to_reduce[-2][index]+=max(last_row[index:index+2])
    del to_reduce[-1]

def find_max_sum(triangle):
    while len(triangle)>1:
        _reduce_triangle(triangle)
    return triangle[0][0]

def _parse_triangle_from_file(data_file):
    triangle=[]
    with open(data_file,'r')as triangle_file:
        for line in triangle_file:
            triangle.append([int(x) for x in line.split()])
    return triangle

triangle=_parse_triangle_from_file('triangle_18.txt')
print find_max_sum(triangle)

