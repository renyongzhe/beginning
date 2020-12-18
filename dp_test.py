def paths(m,n):
    res = [[1]*n]*m
    for i in range(1,m):
        for j in range(1,n):
            res[i][j] = res[i-1][j] + res[i][j-1]
    return res[-1][-1]

paths(10,6)

def path2(m,n):
    if m<= 0 or n<=0 :
        return 0
    dp = [[0]*n]*m
    for i in range(m):
        dp[i][0] = 1
    for j in range(n):
        dp[0][j] = 1
    for i in range(1,m):
        for j in range(1,n):
            dp[i][j] = dp[i-1][j] + dp[i][j-1]
    return dp[m-1][n-1]

path2(12,5)

def jump(n):
    res = list(range(n+1))
    for i in range(1,n+1):
        if i <= 2:
            res[i] = i
        else:
            res[i] = res[i-2] + res[i-1]
    return res[-1]

jump(3)

def fib(n):
        res = list(range(n+1))
        for i in range(n+1):
            if i < 2:
                res[i] = i
            else:
                res[i] = res[i-1] + res[i-2]
        return res[-1]
fib(6)

arr = [-2,1,-3,4,-1,2,1,-5,4]
def maxSubarr(arr):
    target = [ 0 for i in arr]
    maxV = arr[0]
    target[0] = arr[0]
    for i in range(1,len(arr)):
        target[i] = arr[i] + max(target[i-1],0)
        maxV = max(maxV,target[i])
    return maxV,target

aa,cc = maxSubarr(arr)
print(aa)
print(cc)
