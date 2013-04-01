"""
Very simple module to solve linear equation systems in finite fields.


"""
import numpy as np
from gmpy import mpz, invert
import pdb  # XXX: Delete when done
from math import floor

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

def indep_solutions(solved, base):
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

def rat_poly_sing(mat, indep, d1, d2, base, coeff=False):
    ## XXX: Only for singular solutions for now
    coefficients = [0] * (d1 + d2)
    try:
        coefficients[list(indep).index(1)] = mpz(1)
    except ValueError:
        pass
    for (i, v) in enumerate(indep[:-1]):
        if v == 0:
            coefficients[i] = (mat[i][-1] - mat[i][list(indep).index(1)]) % base
    coeff_n = list(enumerate(coefficients[d1-1::-1]))
    coeff_d = list(enumerate(coefficients[:d1-1:-1]))
    if coeff:
        return coeff_n, coeff_d
    f = lambda x: (pow(x, d1, base) + sum([pow(x, a, base) * b for (a,b) in coeff_n])) % base
    g = lambda x: (pow(x, d2, base) + sum([pow(x, a, base) * b for (a,b) in coeff_d])) % base
    return f,g

def evaluate(hostset, points, base):
    values = []
    char_coef = np.poly(hostset)
    for point in points:
        evaluated = np.polyval(char_coef, point) % base
        values.append(mpz(int(evaluated)))
    return values
    
def divide(set1, set2, base):
    values = []
    for (v1, v2) in zip(set1, set2):
        values.append(v1 * v2.invert(base) % base)
    return values

def poly_bounds(set1, set2, m):
    delta = len(set1) - len(set2)
    d1 = int(floor((m+delta)/2.0))
    d2 = int(floor((m-delta)/2.0))
    return d1,d2

def test(fn1, fn2, base):
    ## _very_ naive. still working on root finding
    sol1 = []
    sol2 = []
    for i in range(base):
        if fn1(i) == 0:
            sol1.append(i)
    for i in range(base):
        if fn2(i) == 0:
            if i in sol1:
                sol1.remove(i)
            else:
                sol2.append(i)
    return sol1, sol2

def test_data():
    evpoints = [-1, -2, -3, -4, -5]
    values = [75, 74, 17, 1, 35]
    d1 = 2
    d2 = 3
    base = 97
    return evpoints, values, d1, d2, base

