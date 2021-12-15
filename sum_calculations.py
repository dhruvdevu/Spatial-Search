import numpy as np
import time
from multiprocessing import Pool, freeze_support
from itertools import repeat
from sklearn.utils.extmath import cartesian

d=5
n=20
wVec = np.full(d,n/2)

def term(w, k):
    num = 1
    denom = d
    for i in range(len(w)):
        num *= (np.cos(np.pi*k[i]*w[i]/(2*n)))**2
        denom -= np.cos(np.pi*k[i]/n)
    if denom == 0:
        print(w,k)
    return 0.5*num/denom

def oldterm(w,k):
    denom = d
    for i in range(len(w)):
        denom -= np.cos(np.pi*k[i]/n)
    if denom == 0:
        print(k)
    return 0.5*1/denom

def ktuples():
    karr = np.arange(0,n)
    dupl = np.tile(karr, (d,1))
    cart = cartesian(dupl)
    # remove [00000]
    return cart[1:]
def term_wrapper(k):
    return term(wVec, k)
def oldterm_wrapper(k):
    return oldterm(wVec, k)
def main():
    wVec = np.full(d,n/2)

    tuple_time_start = time.time()
    kList = ktuples()
    tuple_time_end = time.time()

    # multiprocessing
    pool = Pool()
    print("Testing threading:")
    vtime = time.time()
    vals = np.array(pool.starmap(term, zip(repeat(wVec), kList)))
    #vals2 = np.array(pool.map(term_wrapper, kList))
    print("vtime = ", time.time()-vtime)
    oldvtime = time.time()
    oldvals = np.array(pool.starmap(oldterm, zip(repeat(wVec), kList)))
    print("Old vtime = ", time.time()-oldvtime)


    # vals_time_start = time.time()
    # #vals = np.array([term(wVec, k) for k in kList])
    # vals_time_end = time.time()
    # oldvals_time_start = time.time()
    # oldvals = np.array([oldterm(wVec, k) for k in kList])
    # oldvals_time_end = time.time()
    sum_time_start = time.time()
    sum = ((2/n) ** d)*np.sum(vals)
    sum_time_end = time.time()
    oldsum_time_start = time.time()
    oldsum = ((1/n) ** d)*np.sum(oldvals)
    oldsum_time_end = time.time()
    print(sum, oldsum)
    print("Tuple time = ", tuple_time_end-tuple_time_start)
    # print("vals time = ", vals_time_end-vals_time_start)
    # print("oldvals time = ", oldvals_time_end-oldvals_time_start)
    print("sum time = ", sum_time_end-sum_time_start)
    print("old sum time = ", oldsum_time_end-oldsum_time_start)

if __name__ == "__main__":
    main()
