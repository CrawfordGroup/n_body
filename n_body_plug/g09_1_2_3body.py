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
body_list += [1]
# two-body
body_list += [2]
# three-body
body_list += [3]
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

method_list += ['b3lyp']
#method_list += ['ccsd']
#method_list += ['cc2']

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
#       Example: db['cc2'][2]['job_status'] = OrderedDict([('1','not_started'),('2','not_started')])
        body = n_body_dir(n)
#       Example: n_body_dir(2) = 'two_body'
        while job_list:
            job,val = job_list.popitem(last=True)
#           Example: job = 2 and val = not_started
#           Actually pops the entry out of the job_list so that while job_list actually ends
#            if val == 'not_started': # For now, only start not_started jobs. Start dead jobs later
            if (val == 'not_started') or (val == 'error'):
                print(job)
                print('{}/{}/{}'.format(method,body,job))
                os.chdir('{}/{}/{}'.format(method,body,job))
#               Example: os.chdir('{}/{}/{}'.format('cc2','two_body',2)
                subprocess.call(['g09', 'input.com', 'input.log'])
                os.chdir('../../..')
        db.close()

