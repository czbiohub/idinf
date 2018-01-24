#!/usr/bin/env python3

import sys
import networkx
import numpy as np
import scipy.sparse as sp
from cvxpy import *
from collections import defaultdict

def notes_eg():
    ''' In this example, we have three species with multimapping reads.
     | n   | species 1 | species 2 | species 3 |
     |-----|-----------|-----------|-----------|
     | 100 |     X     |           |           |
     | 10  |           |     X     |           |
     | 10  |           |           |     X     |
     | 40  |     X     |     X     |           |
     | 30  |     X     |     X     |     X     |
    '''
    D = sp.diags([100, 10, 10, 40, 30])
    X_dense = [[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1],
               [1, 1, 0],
               [1, 1, 1]]
    X = sp.csr_matrix(X_dense)
    q_exact_optimum = np.array([150/176,15/176,11/176])
    return X, D, q_exact_optimum


def two_species_eg():
    ''' In this example, we have two species with multimapping reads,
     only one of which is present.
     | n   | species 1 | species 2 |
     |-----|-----------|-----------|
     | 1   |     1     |     1     |
     | 1   |     1     |     0     |
    '''
    D = sp.diags([1, 1])
    X_dense = [[1, 1],
               [1, 0]]
    X = sp.csr_matrix(X_dense)
    q_exact_optimum = np.array([1,0])
    return X, D, q_exact_optimum


def test_eg(X, D, q_exact_optimum):

    n_read_types, n_species = X.shape
    print(f"Testing example with {n_species} species, {n_read_types} read types and {D.sum()} reads.")

    q_opt = cvx_solve(X, D)

    if (np.allclose(q_opt, q_exact_optimum, atol=1e-6)):
        print("Test passes.")
    else:
        print("Convex solver did not find the solution.")

def cvx_solve(X, D, verbose=False):
    '''
    :param X: Matrix whose rows are patterns of mappings and columns are species
    :param D: Diagonal matrix recording counts of each pattern of mapping
    :param verbose: Verbosity for the convex solver.
    :return: Array of optimal probabilities.
    '''
    n = X.shape[1]
    d = diag(D)

    q = Variable(n)
    Y = X * q

    objective_log = 0
    for i in range(X.shape[0]):
        row = log(Y[i]) * d[i]
        objective_log += row

    prob = Problem(Maximize(objective_log), [q >= 0, sum(q) == 1])
    prob.solve(verbose=verbose)

    q_array = np.array(q.value).ravel()

    return q_array


def read_proportion(X, D):
    normalizer = sp.diags((1 / X.sum(axis=1)).A1)
    read_proportion = normalizer.dot(D.dot(X)).sum(axis=0).A1
    read_proportion = read_proportion / read_proportion.sum()
    return read_proportion


def grad(X, D, q):
    return X.T.dot(D.dot(1 / X.dot(q)))


def proj_grad(X, D, q):
    g = grad(X, D, q)
    return g - g.mean()


def proj_grad_descent(alpha, q_0, proj_grad_fun, max_iter):
    iterates = np.zeros((max_iter, len(q_0)))
    q = q_0
    for i in range(max_iter):
        iterates[i, :] = q
        x = q + alpha * proj_grad_fun(q)
        x = x * (x > 0)
        q = x / x.sum()
    return q, iterates


def test_eg_proj_gradient(X, D, q_exact_optimum):
    n_read_types, n_species = X.shape
    alpha = 1e-1 / D.sum()
    max_iter = 1000

    q_0 = read_proportion(X, D)

    print(f"Testing example with {n_species} species, {n_read_types} read types and {D.sum()} reads.")
    q_opt, iterates = proj_grad_descent(alpha, q_0, lambda x: proj_grad(X, D, x), max_iter)
    if (np.allclose(q_opt, q_exact_optimum)):
        print("Test passes.")
    else:
        print("Grad descent did not converge to the optimum in {max_iter} steps.")






