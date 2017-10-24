#!/usr/bin/env python

import shutil
import shelve
import os
import subprocess

## Submission script
#  This script easily allows you to submit subsets of the n_body jobs broken
#  down over approximation levels
#  Jobs are submitted in reverse order by default because we're assuming the 
#  solute(s) were the first fragments

## List of jobs to be submitted
#  Starts out empty
body_list = []
method_list = []

# one-body
#body_list += [1]
# two-body
body_list += [2]
# three-body
#body_list += [3]
# four-body
#body_list += [4]
# five-body
#body_list += [5]
# six-body
#body_list += [6]
# seven-body
#body_list += [7]
# eight-body
#body_list += [8]
# nine-body
#body_list += [9]
# ten-body
#body_list += [10]
# eleven-body
#body_list += [11]

#method_list += ['b3lyp']
#method_list += ['ccsd']
method_list += ['cc2']

## Submit file for jobs containing the solute molecule
#  Default name is solute
solute = 'solute'
## Submit file for jobs containing only solvent molecules
#  Default name is solvent
solvent = 'solvent'

## Solute fragments
#  Default is fragment 1
#  Everything else is assumed to be solvent
solute_frags = set([1])

def n_body_dir(n):
    num_2_word = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven',
                  'eight', 'nine', 'ten', 'eleven', 'twelve', 'thirteen', 'fourteen']

    return('{}_body'.format(num_2_word[n]))


for method in method_list:
    for n in body_list:
        db = shelve.open('database')
        job_list = db[method][n]['job_status']
        body = n_body_dir(n)
        while job_list:
            job,val = job_list.popitem(last=True)
            print(job)
            os.chdir('{}/{}/{}'.format(method,body,job))
            subprocess.call(['psi4', 'input.dat', 'output.dat'])
            os.chdir('../../..')

