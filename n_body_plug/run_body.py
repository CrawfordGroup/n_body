#!/usr/bin/env python
#import shutil
import shelve
import os
import subprocess
import numpy as np

# grab input args at command line
if __name__ == "__main__":
    import sys
    method = sys.argv[1] # level of theory, eg cam-b3lyp (one at a time)
    n = int(sys.argv[2]) # body to submit, eg 4 (one at a time, may want to extend)
    if len(sys.argv) > 3: # if I give args for splitting the job_list into multiple, get them
        n_split = int(sys.argv[3]) # split job_list into n_split lists
        sub = int(sys.argv[4]) # run the sub-th list from above split
    else:
        n_split = 1 # submit all n-body jobs at once
        sub = 1 # ditto above

def n_body_dir(n):
    '''Returns string of spelled-out number for directory crawling'''
    num_2_word = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven',
                  'eight', 'nine', 'ten', 'eleven', 'twelve', 'thirteen', 'fourteen']
    return('{}_body'.format(num_2_word[n]))


db = shelve.open('database')
job_list = db[method][n]['job_status'] # total job list
for job,val in list(job_list.items()): # iterate over list.items so size doesn't change
    if val == 'complete': # only submit running or error'd jobs
        del job_list[job]
job_list = list(job_list.keys()) # list of all unsubbed jobs
j_lists = np.array_split(job_list, n_split) # split into equal(ish) np.arrays
use = j_lists[sub - 1] # jobs to run: only the sub-th array
body = n_body_dir(n) # body name for directory crawling
while (len(use) > 0): # run all in use list in reverse order (solute probably first)
    job, use = use[-1], use[:-1] # pull the job out of the array to execute
    print(job)
    print('{}/{}/{}'.format(method,body,job))
    os.chdir('{}/{}/{}'.format(method,body,job))
#   Example: os.chdir('{}/{}/{}'.format('cc2','two_body',2)
    subprocess.call(['g09', 'input.com', 'input.log'])
    os.chdir('../../..')
db.close()

