def is_palin(string):
    start, end =0, len(string)-1
    while end > start:
        if string[start] != string[end]:
                return False
        start += 1
        end -= 1
    return True

def palindrome():
    num1=999
    result_arr = []
    while num1 > 100:
        num2=990
        if num2>num1:
            num2=num1-(num1&11)
        while num2>109:
            if is_palin(str(num1 * num2)):
                result_arr.append(num1 * num2)
            num2 -= 11
        num1 -= 1
    result_arr.sort()
    print result_arr[len(result_arr)-1]

palindrome()
