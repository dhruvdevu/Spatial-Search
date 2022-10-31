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


# In[ ]:


psums = []
apsums = []
for i in [47]:
    wVec = np.full(d,n/2)
    n=i
    d=5
    ps = eval_periodic_sum()
    s = eval_aperiodic_sum(wVec)
    print(n, ps, s)
    apsums += [s]
    psums += [ps]

print(psums)
print(apsums)
plt.plot(psums)
plt.plot(apsums)
plt.show()
