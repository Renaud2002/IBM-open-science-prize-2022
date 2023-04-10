
# Author: Edison Murairi
# Date: February 21st, 2023
# Given a stabilizer S, find a Pauli opertor in N(S) with miniamal weight d

#%%
import diag_helpers as helpers
import itertools as it
from diag_optimization import weight
import numpy as np
import galois
GF = galois.GF(2)


from tableau_operations import makeTableauMatrix

def generators_strings(x,z):

    return helpers.tableau_to_pstrings(x,z)

def template(k,n,r):

    return it.combinations(range(r), r=k)

def apply_template(nullspace, template):

    res = GF([0]*nullspace.shape[1])

    for t in template:

        res += nullspace[t]

    return res

def brute_force_optimize(nullspace):

    r = nullspace.shape[0]
    n = nullspace.shape[1]//2
    k = 1

    curr_vector = nullspace[0]
    curr_weight = weight(curr_vector)

    while k <= 3:

        templts = template(k, n, r)

        for templ in templts:

            tmp_vector = apply_template(nullspace, templ)
            tmp_weight = weight(tmp_vector)

            if tmp_weight < curr_weight:

                curr_vector = tmp_vector
                curr_weight = tmp_weight

        k += 1

        # print("n = {0}, r = {1}, curr_weight = ".format(n, r), curr_weight)

        # if curr_weight <= r/2 - 4:
        #     return curr_vector
    
    # print(f"r = {n}, weight = {curr_weight}")

    return curr_vector

#%% My second implementation of Brute force. This time, I search through the whole normalizer
def commute(p1, p2):

    assert len(p1) == len(p2)
    n = len(p1)//2

    return int(p1[:n].dot(p2[n:]) + p1[n:].dot(p2[:n])) == 0

def commute_with_others(p, others):
    for other in others:
        if not commute(p, other):
            return False

    return True
# %%
def normalizer(T):

    n = T.shape[1]//2
    candidates = it.product([0,1], repeat = 2*n)

    for candidate in candidates:

        candidate = GF(candidate)

        if commute_with_others(candidate, T):
            yield candidate

# %%
def search_normalizer(x,z):
    T = makeTableauMatrix(z,x)
    normal = normalizer(T)
    #normal = list(normal)

    currp = T[0]
    currw = weight(currp)


    for nextp in normal:

        nextw = weight(nextp)
        if nextw < currw and nextw !=0:
            currp = nextp
            currw = nextw
            

    # print("weight: ", currw)

    return currp
