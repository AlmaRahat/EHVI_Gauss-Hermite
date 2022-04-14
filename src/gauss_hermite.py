#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions to generate Gauss-Hermite nodes for single and multivariate Gaussian 
Distirbutions. 

This code is an adaptation of the following blog post:
https://www.r-bloggers.com/2015/09/notes-on-multivariate-gaussian-quadrature-with-r-code/

:Author:
   Alma Rahat   <a.a.m.rahat@swansea.ac.uk>
:Date:
   14 April 2022
:Copyright:
   Copyright (c)  Alma Rahat, Swansea University, 2022
"""

# imports
import numpy as np

def hermite(points, z):
    p1 = 1/np.pi**0.4
    p2 = 0
    for i in range(1, points+1):
        p3 = p2
        p2 = p1
        p1 = z * np.sqrt(2/i) * p2 - np.sqrt((i - 1)/i) * p3
    pp = np.sqrt(2 * points) * p2
    return np.array([p1, pp])


def gauss_hermite(points, interlim=50):
    x = np.zeros(points)
    w = np.zeros(points)
    m = np.int(np.floor((points + 1)/2))
    for i in range(1, m+1):
        if i == 1:
            z = np.sqrt(2 * points +1) - 2 * (2 * points + 1)**(-1/6)
        elif i == 2:
            z = z - np.sqrt(points)/z
        elif (i == 3 or i == 4):
            z = 1.9 * z - 0.9 *x[i-1-2]
        else:
            z = 2 * z - x[i-1-2]
        # Newton-Raphson loop
        for j in range(interlim):
            z1 = z
            p = hermite(points, z)
            z = z1 - p[0]/p[1]
            if np.abs(z1 - z)<1e-15:
                break
        if j == interlim - 1:
            warnings.warn("iteration limit reached!")
        x[i-1] = z
        x[points - i] = -z
        f = 2/p[1]**2
        w[i-1] = f
        w[points - i] = f
    return x * np.sqrt(2), w/np.sum(w)

def mGauss_hermite(n, mu, sigma, prune=None):
    dm = len(mu)
    gh = gauss_hermite(n)
    l1 = [gh[0]]*dm 
    l2 = [gh[1]]*dm 
    x = np.array(np.meshgrid(*l1)).transpose().reshape((-1, dm))
    w = np.prod(np.array(np.meshgrid(*l2)).transpose().reshape(-1,dm),
                    axis=1)
    if prune is not None:
        qwt = np.quantile(w, prune)
        inds = np.where(w>qwt)[0]
        x = x[inds]
        w = w[inds]
    eigval, eigvec = np.linalg.eig(sigma)
    rot = np.dot(eigvec, np.diag(np.sqrt(eigval)))
    x = np.dot(rot,x.transpose()).transpose()
    return x+mu, w



