#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python wrapper the analytical EHVI for 2D and 3D using the code available in:
https://liacs.leidenuniv.nl/~csmoda/index.php?page=code

:Author:
   Alma Rahat   <a.a.m.rahat@swansea.ac.uk>
:Date:
   14 April 2022
:Copyright:
   Copyright (c)  Alma Rahat, Swansea University, 2022
"""

# imports
import numpy as np
import subprocess as cmd

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

def line_appender(filename, line):
    with open(filename, 'a') as f:
        f.writelines(line)
        
def generate_txt_input(pf_filename, r, mu, std, filename="test.txt", multiplier=-1):
    prob = pd.read_csv(pf_filename, header=None, delimiter=" ", dtype=float)
    pf= prob.to_numpy()
    leng = pf.shape[0]
    np.savetxt(filename, multiplier*pf, delimiter=" ", fmt="%.15f")
    line_prepender(filename, str(leng)+"\n")
    rapp = "\n" + np.array2string(multiplier*r, precision=15).strip('[]')
    line_appender(filename, rapp)
    muapp = "\n" +np.array2string(multiplier*mu, precision=15).strip('[]')
    line_appender(filename, muapp)
    stdapp = " " +np.array2string(std, precision=15).strip('[]')
    line_appender(filename, stdapp)
    
def generate_txt_input_pf(pf, r, mu, std, filename="test.txt", multiplier=-1):
    # multiplier = -1 for minimisation
    leng = pf.shape[0]
    np.savetxt(filename, multiplier*pf, delimiter=" ", fmt="%.15f")
    line_prepender(filename, str(leng)+"\n")
    rapp = "\n" + np.array2string(multiplier*r, precision=15).strip('[]')
    line_appender(filename, rapp)
    muapp = "\n" +np.array2string(multiplier*mu, precision=15).strip('[]')
    line_appender(filename, muapp)
    stdapp = " " +np.array2string(std, precision=15).strip('[]')
    line_appender(filename, stdapp)

def run_ehvi(code_base, conf_filename):
    a = cmd.check_output(['cp', conf_filename, code_base])
    b = cmd.check_output([code_base+"EHVI", conf_filename])
    return float(b)