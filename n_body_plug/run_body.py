#!/usr/bin/env python
#import shutil
import shelve
import os
import subprocess
import numpy as np

# List of jobs to be submitted
if __name__ == "__main__":
    import sys
    method = sys.argv[1]
    n = int(sys.argv[2])
    if len(sys.argv) > 3:
        n_split = int(sys.argv[3])
        sub = int(sys.argv[4])
    else:
        n_split = 1
        sub = 1

#body_list = [7]
#method_list = ['cam-b3lyp']
#n_split = 3 # split total jobs list into n_split-many arrays of jobs
#sub = 1 # submit the sub-th array of the n_split arrays created

def n_body_dir(n):
    num_2_word = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven',
                  'eight', 'nine', 'ten', 'eleven', 'twelve', 'thirteen', 'fourteen']
    return('{}_body'.format(num_2_word[n]))


db = shelve.open('database')
job_list = db[method][n]['job_status']
for job,val in list(job_list.items()): # iterate over list.items so size doesn't change
    if val == 'complete': # only submit running or error'd jobs
        del job_list[job]
job_list = list(job_list.keys()) # get a list of all unsubbed jobs
j_lists = np.array_split(job_list, n_split) # split into equal(ish) np.arrays
use = j_lists[sub - 1] # we're only using the sub-th array
body = n_body_dir(n)
while (len(use) > 0):
    job, use = use[-1], use[:-1] # pull the job out of the array to execute
    print(job)
#    print('{}/{}/{}'.format(method,body,job))
    os.chdir('{}/{}/{}'.format(method,body,job))
#   Example: os.chdir('{}/{}/{}'.format('cc2','two_body',2)
#    subprocess.call(['g09', 'input.com', 'input.log'])
    os.chdir('../../..')
db.close()

