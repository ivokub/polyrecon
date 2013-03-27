"""
Very simple module to solve linear equation systems in finite fields.


"""
import numpy as np
from gmpy import mpz, invert

def solve_bottom(mat, base):
    localmat = mat.copy() 
    n, m = mat.shape
    for i in range(n-1):
        localmat[i] = inverse_array(localmat[i], base, i)
        for j in range(i+1,n):
            localmat[j] += -localmat[j][i] * localmat[i]
    return mod_matrix(localmat, base)

def solve_top(mat, base):
    localmat = mat.copy()
    n, m = mat.shape
    for i in range(1,n)[::-1]:
        localmat[i] = inverse_array(localmat[i], base, i)
        for j in range(i)[::-1]:
            localmat[j] += -localmat[j][i] * localmat[i]
    return mod_matrix(localmat, base)

def solve(mat, base):
    localmat = mat.copy()
    localmat = solve_bottom(localmat, base)
    localmat = solve_top(localmat, base)
    return localmat
    
def mod_array(arr, base):
    return arr % base

def mod_matrix(mat, base):
    return mat % base

def inverse_array(arr, base, pos = 0):
    inverse = invert(arr[pos], base)
    return mod_array(arr * inverse, base)

def create_array(arr):
    mpzarr = [mpz(i) for i in arr]
    a = np.array(mpzarr)
    return a
    
def create_mat(mat):
    npmat = np.array([create_array(vec) for vec in mat])
    return npmat

def create_equations(evpoints, values, d1, d2, base):
    mat = []
    for (ev, val) in zip(evpoints, values):
        vec = []
        for i in range(d1 - 1, -1, -1):
            vec.append(pow(ev, i, base))
        for i in range(d2 - 1, -1, -1):
            vec.append(-val * pow(ev, i, base))
        vec.append(val * pow(ev, d2, base) - pow(ev,d1, base))
        mat.append(vec)
    return mod_matrix(create_mat(mat), base)

def rat_poly(solved, base):
    fil = lambda x: x != 0
    ## TODO: Consider if solution is not sigular and how to compute all
    ## solutions
    independent = np.array([0] * solved.shape[1])
    dependent = np.array([0] * solved.shape[1])
    for row in solved:
        dependent |= map(fil, row) & independent
        independent |= map(fil, row)
    dependent[-1] = 0    # answer may be dependent, we dont care
    return dependent, independent

def test_data():
    evpoints = [-1, -2, -3, -4, -5]
    values = [75, 74, 17, 1, 35]
    d1 = 2
    d2 = 3
    base = 97
    return evpoints, values, d1, d2, base

