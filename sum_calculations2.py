#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import time
import matplotlib.pyplot as plt
from multiprocessing import Pool, freeze_support
from itertools import repeat
from sklearn.utils.extmath import cartesian


# In[2]:


d=5
n=23
wVec = np.full(d,n/2)


# In[3]:


def term(w, k):
    num = 1
    denom = d
    for i in range(len(w)):
        num *= (np.cos(np.pi*k[i]*(2*w[i]-1)/(2*n)))**2
        denom -= np.cos(np.pi*k[i]/n)
    if denom == 0:
        print(w,k)
    return 0.5*num/denom


# In[4]:


def fullterm(w, k):
    num = 1
    denom = d
    for i in range(len(w)):
        num *= (np.cos(np.pi*k[i]*(2*w[i]-1)/(2*n)))**2
        denom -= np.cos(np.pi*k[i]/n)
    if denom == 0:
        print(w,k)
    return (num, denom, 0.5*num/denom)


# In[5]:


def oldterm(k):
    denom = d
    for i in range(d):
        denom -= np.cos(np.pi*k[i]/n)
    if denom == 0:
        print(k)
    return 0.5*1/denom


# In[6]:


def ktuples(r):
    karr = np.arange(0,r)
    dupl = np.tile(karr, (d,1))
    cart = cartesian(dupl)
    # remove [00000]
    return cart[1:]


# In[7]:


def eval_periodic_sum():
    kList = ktuples(n)
    pool = Pool()
    pvals = np.array(pool.map(oldterm, kList))
    periodic_sum = ((1/n) ** d)*np.sum(pvals)
    return periodic_sum


# In[8]:


def eval_aperiodic_sum(w):
    kList = ktuples(n)
    pool = Pool()
    vals = np.array(pool.starmap(term, zip(repeat(w), kList)))
    vsum = ((2/n) ** d)*np.sum(vals)
    return vsum


# In[10]:


n=21
d=5
eval_periodic_sum()


# In[10]:


n=20
wVec = np.full(d,n/2)
start_time = time.time()
s = eval_aperiodic_sum(wVec)
print(s)
print(time.time()-start_time)


# In[11]:


print((25*((5/2)**5))/60)
print(np.full(5,int(21/2)))


# In[ ]:


psums = []
apsums = []
for i in [49,51]:
    wVec = np.full(d,n/2)
    n=i
    d=5
    ps = eval_periodic_sum()
    s = eval_aperiodic_sum(wVec)
    print(n, ps, s)
    apsums += [s]
    psums += [ps]


# In[11]:


psums = [0.1246195030505344, 0.12415854206923803, 0.12336459907023252, 0.12270531347837141, 0.12146128093489232, 0.1203083486927991, 0.11994270934702529, 0.11983319181216524]
apsums = [0.13501651292491895,0.15856251639028932,0.1322090856592393, 0.13071932948957735, 0.1279704600891266,0.1385064780439989, 0.12462302411104727, 0.1245096927963745]


# In[12]:


plt.plot(psums)
plt.plot(apsums)
plt.show()


# In[21]:


n=20
print(eval_periodic_sum())


# In[22]:


n=40
wVec = np.full(d,n/2)
start_time = time.time()
s = eval_aperiodic_sum(wVec)
print(s)
print(time.time()-start_time)


# In[8]:


wVec = np.full(d,n/2)
kList = ktuples(n)
oneVec = np.full(d, 1)
pool = Pool()
sums = []
for r in range(n):
    w = np.array([r+1,(n+1)/2,(n+1)/2,(n+1)/2,(n+1)/2])
    vtime = time.time()
    vals = np.array(pool.starmap(term, zip(repeat(w), kList)))
    vtime = time.time()-vtime
    vsum = ((2/n) ** d)*np.sum(vals)
    sums += [vsum]
    print("Point ", w, " took time ", vtime, " with aperiodic sum = ", vsum)
sums = np.array(sums)

oldvtime = time.time()
oldvals = np.array(pool.map(oldterm, kList))
oldsum = ((1/n) ** d)*np.sum(oldvals)
print("Periodic sum = ", oldsum, ", took time ", time.time()-oldvtime)
plt.plot(sums, color='tab:blue', label="Aperiodic sums")
#plt.plot(oldsum*np.ones(1), color='tab:orange', label="Periodic sum")
plt.ylabel("$S_{1,5}$")
plt.xlabel("x")
plt.title("Sum variation with position along (x,1,1,1,1). n = "+str(n))
plt.legend()
plt.show()


# In[23]:


start_time = time.time()
ps = eval_periodic_sum()
print(ps)
print(time.time()-start_time)


# In[13]:


#check for odd
n=21
wVec = np.full(d,(n-1)/2)
start_time = time.time()
s2 = eval_aperiodic_sum(wVec)
print(s2)
print(time.time()-start_time)
start_time = time.time()
ps2 = eval_periodic_sum()
print(ps2)
print(time.time()-start_time)


# In[58]:


def trit_func(i):
    
    if i%3 == 0:
        return 1.0/6.0
    elif i%3 == 1:
        return 3.0/16.0
    else:
        return 1.0/8.0
        
def trit_strings(d):
    digits = np.arange(3)
    dupl = np.tile(digits, (d,1))
    cart = cartesian(dupl)
    return cart

def trit_func_prod(arr):
    prod = 1
    for a in arr:
        prod *= trit_func(a)
    return prod


trit_strings(5)
    


# In[59]:


trit_sum = 0
for digit_string in trit_strings(5):
    trit_sum += trit_func_prod(digit_string)
print(32*trit_sum)
    
    


# In[60]:


1/32


# In[9]:


ktuples(2)


# In[10]:


# See where the majority of terms in the sum are coming from
n=21
wVec = np.full(d,(n+1)/2)
kList = ktuples(n)
oneVec = np.full(d, 1)
pool = Pool()
sums = []

# Compute the indices in klist at which the radius increases
# Never mind, they are interspersed randomly.



for r in range(1,n+1):
    kList = ktuples(r)
    vtime = time.time()
    vals = np.array(pool.starmap(term, zip(repeat(wVec), kList)))
    vtime = time.time()-vtime
    vsum = ((2/n) ** d)*np.sum(vals)
    sums += [vsum]
    print("Sum upto radius ", r-1, " took time ", vtime, " with aperiodic sum = ", vsum)
sums = np.array(sums)

kList = ktuples(n)
oldvtime = time.time()
oldvals = np.array(pool.map(oldterm, kList))
oldsum = ((1/n) ** d)*np.sum(oldvals)
print("Periodic sum = ", oldsum, ", took time ", time.time()-oldvtime)
plt.plot(sums, color='tab:blue', label="Aperiodic sums")
#plt.plot(oldsum*np.ones(1), color='tab:orange', label="Periodic sum")
plt.ylabel("$S_{1,5} upto radius r$")
plt.xlabel("r")
plt.title("Sum upto radius r. n = "+str(n))
plt.legend()
plt.show()


# In[9]:


n=21
wVec = np.full(d,(n+1)/2)
kList = ktuples(n)
oneVec = np.full(d, 1)
pool = Pool()
sums = []

for r in range(1, n+1):
    kList = ktuples(r)
    vtime = time.time()
    vals = np.array(pool.map(oldterm, kList))
    vtime = time.time()-vtime
    vsum = ((1/n) ** d)*np.sum(vals)
    sums += [vsum]
    print("Sum upto radius ", r-1, " took time ", vtime, " with aperiodic sum = ", vsum)
sums = np.array(sums)

kList = ktuples(n)
oldvtime = time.time()
oldvals = np.array(pool.map(oldterm, kList))
oldsum = ((1/n) ** d)*np.sum(oldvals)
print("Periodic sum = ", oldsum, ", took time ", time.time()-oldvtime)
plt.plot(sums, color='tab:blue', label="Periodic sum")
#plt.plot(oldsum*np.ones(1), color='tab:orange', label="Periodic sum")
plt.ylabel("$S_{1,5} upto radius r$")
plt.xlabel("r")
plt.title("Sum upto radius r. n = "+str(n))
plt.legend()
plt.show()


# In[26]:


n=21
wVec = np.full(d,(n+1)/2)
kList = ktuples(n)
oneVec = np.full(d, 1)
terms = [0]
nums = []
denoms = []
for r in range(n):
    k = (1+r)*oneVec
    num, denom, val = fullterm(wVec, k)
    terms += [val]
    nums += [num]
    denoms += [1/denom]
    
terms = np.array(terms)
nums = np.array(nums)
denoms = np.array(denoms)
print(terms)
print(nums)
print(denoms)
plt.plot(terms, color='tab:blue', label="Aperiodic sums")
#plt.plot(oldsum*np.ones(1), color='tab:orange', label="Periodic sum")
plt.ylabel("$Term from the sum$")
plt.xlabel("Position along (1,1,1,1,1)")
plt.title("Term variation , with w at center and at k = x*(1,1,1,1,1). n = "+str(n))
plt.legend()
plt.show()

plt.plot(nums, color='tab:green', label="Numerators")
#plt.plot(oldsum*np.ones(1), color='tab:orange', label="Periodic sum")
plt.ylabel("Numerator ")
plt.xlabel("Position along (1,1,1,1,1)")
plt.title("Numerator variation with from the sum, with w at center and at k = x*(1,1,1,1,1). n = "+str(n))
plt.legend()
plt.show()

plt.plot(denoms, color='tab:green', label="Denominators")
#plt.plot(oldsum*np.ones(1), color='tab:orange', label="Periodic sum")
plt.ylabel("1/Denominator ")
plt.xlabel("Position along (1,1,1,1,1)")
plt.title("1/denominator variation with from the sum, with w at center and at k = x*(1,1,1,1,1). n = "+str(n))
plt.legend()
plt.show()


# In[35]:


unperturbed = []
perturbed = []
delta = 19
n=41
s1 = 0
s2 = 0
for k in range(n):
    denom = 2-np.cos(np.pi*k/n)
    val1 = (np.cos(np.pi*k/2) ** 2)/denom
    val2 = (np.cos(np.pi*k/2 + np.pi*k*delta/n)**2)/denom
    unperturbed += [val1]
    perturbed += [val2]
    s1 += val1
    s2 += val2

print(s1)
print(s2)
plt.plot(unperturbed, color='tab:green', label="Unperturbed")
plt.plot(perturbed, color='tab:blue', label="Perturbed by 1")
plt.xlabel("k")
plt.ylabel("$\cos^2 \pi k (2w-1)/2n$")
#plt.legend()
plt.show()

    


# In[23]:


unperturbed = []
perturbed = []
delta = 20
n=41
s1 = 0
s2 = 0
for k in range(n):
    denom = 2-np.cos(np.pi*k/n)
    val1 = (np.cos(np.pi*k/5) ** 2)/denom
    val2 = (np.cos(np.pi*k/5 + np.pi*k*delta/n)**2)/denom
    unperturbed += [val1]
    perturbed += [val2]
    s1 += val1
    s2 += val2

print(s1)
print(s2)
plt.plot(unperturbed, color='tab:green', label="Unperturbed")
plt.plot(perturbed, color='tab:blue', label="Perturbed by 1")
plt.xlabel("k")
plt.ylabel("$\cos^2 \pi k (2w-1)/2n$")
#plt.legend()
plt.show()


# In[15]:


oldterms = [oldterm(i*np.ones(5)) for i in range(1,n)]
plt.plot(oldterms)
2*oldterm(0.7*n*np.ones(5))


# In[35]:


print(0.124/(21 ** 5 - 1))
print(oldterms)


# In[32]:


len(kList)


# 20 ** 5
# 

# In[34]:


21 ** 5 -1


# In[40]:


print(cartesian([np.arange(0,n,2),np.arange(0,n),np.arange(0,n),np.arange(0,n),np.arange(0,n)]))


# In[42]:


#Compare odd and even sums.
kList = ktuples(n)
oneVec = np.full(d, 1)
pool = Pool()
sums = []
evenkList = cartesian([np.arange(0,n,2),np.arange(0,n),np.arange(0,n),np.arange(0,n),np.arange(0,n)])[1:]
oddkList = cartesian([np.arange(1,n,2),np.arange(0,n),np.arange(0,n),np.arange(0,n),np.arange(0,n)])

oldvtime = time.time()
oldvals = np.array(pool.map(oldterm, kList))
oldEvenVals = np.array(pool.map(oldterm, evenkList))
oldOddVals = np.array(pool.map(oldterm, oddkList))
oldsum = ((1/n) ** d)*np.sum(oldvals)
oldEvenSum = ((1/n) ** d)*np.sum(oldEvenVals)
oldOddSum = ((1/n) ** d)*np.sum(oldOddVals)

print("Periodic sum = ", oldsum, ", took time ", time.time()-oldvtime)
print("Periodic even sum = ", oldEvenSum)
print("Periodic odd sum = ", oldOddSum)


# In[ ]:




