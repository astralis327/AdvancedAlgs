import math
import random


def prime_test(N, k):
    # This is main function, that is connected to the Test button. You don't need to touch it.
    mod_exp(2,4,6)
    return fermat(N, k), miller_rabin(N, k)


def mod_exp(x, y, N):
    if y == 0: return 1                     #base case
    z = mod_exp(x, math.floor(y/2), N)
    if y % 2 == 0: return (z**2) % N
    else: return (x*(z**2)) % N             #Complexity is O(n^3)



def euclid_extended(a,b):
    if b == 0: return (1,0,a)
    xprime, yprime, d = euclid_extended(b, a % b)
    return (yprime, xprime - math.floor(a/b) * yprime, d) #Complexity is O(n^3)

def fprobability(k):
    return 1 - pow(1 / 2, k)


def mprobability(k):
    # You will need to implement this function and change the return value.   
    return 1-pow(1/2, k)


def fermat(N, k):
    #Direct application of Fermat's Little Theorem
    for trial in range(k):
        rand_int = random.randint(1, N)     #Generate the test base number
        if mod_exp(rand_int, N-1, N) != 1:  #If a^(N-1) mod N does not equal 1, the test fails
            return 'composite'
    return 'Probably prime'             #Complexity is k*(n^3+1+1)= O(n^4).


def miller_rabin(N, k):

    # Miller Rabin Alg steps:
    # 1) find n-1 = 2^i * m
    # 2) choose a: 1 < a < n-1
    # 3) compute b_0 = a^m (mod n), b_i = (b_(i-1))^2

    #1+1+ log_2(N-1)(2) + k*
    #Step 1
    i = 0               #1
    m = N-1             #1
    while m % 2 == 0:       #This is finding the lowest number a can be raised to and still return 1 when mod N
        m = m >> 1          #log2(N-1)(1+1) Steps
        i = i + 1   #number of times we can take the square root of a^(N-1)

    for iteration in range(k):  #k steps
        maybePrime = False      #1
        #Step 2
        a = random.randint(2, N - 2)    #1
        #Step 3
        x = mod_exp(a, m, N) #a^m (mod N)   #n^3
        if x == 1 or x == N - 1:            #1
            maybePrime = True               #1
            continue        #find another base to test, 1 step
        for j in range(i - 1):              #i-1 steps
            x = mod_exp(x, 2, N)    #Else, keep squaring a^m until a^m (mod N) = N-1. n^3 steps
            if x == N - 1:          #1
                maybePrime = True   #1      1+1+2log2(N-1)+k(1+1+n^3+1+1+1+(i-1)(n^3 +1+1+1) + 1)+1+
                break               #1
        if maybePrime == False: return 'composite' #1
    return 'Probably prime' #1
#
 #1+1+2log2(N-1)+k(1+1+n^3+1+1+1+(i-1)(n^3 +1+1+1) + 1)+1
 #k and i can be simplified to n, so the largest term in the sum above is n^5. 






















