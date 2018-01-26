##########################
### DATABASE STRUCTURE ###
##########################
#database = {
#    'initialized'     : bool, #True after database initialization
#    'inputs_generated': bool, #True after generation of all inputs
#    'jobs_complete'   : bool, #True after all jobs have run
#    'results_computed': bool, #True after final output is produced
#    'num_threads      : int,  #Number of threads for parallel jobs (currently only for g09 interface)
#    'options' = { #Holds user specified options
#        'distributed' : bool, #Generate directories if True, else run in single executable
#        'n_body_func' : function, #Type of job requested (i.e. property)
#        'methods' : list, #Holds strings of all methods requested
#    }
#    method = { #One for each method requested (i.e. 'cc2')
#        'n_body_max': int #User set upper bound for n_body jobs
#        'results': list, #desired computed quantities
#        n_body_level = { #One of each n <= n_body_max (an int, i.e. 1)
                          #Also ghost jobs (a str, i.e. '2-1r' [2 body, 1 real])
#            'num_jobs_complete': int, #
#            'total_num_jobs'   : int, #total number of jobs for n
#            result = {} #dict containing information for each quantity of interest
#            result['raw_data'] = 
#            result['cooked_data'] = 
#            result['correction'] = 
#            result[result] = #contains the total quantity at n_body_level
#            'job_status' = { #OrderedDict of current job status
#                # job_id # : # status #, #Specific examples of possible values follow
#                  '1'      : 'not_started', #No output file (output.dat) present
#                  '1_2'    : 'running', #No 'PSI4 exiting successfully' in output
#                  '1_2_3'  : 'complete' #'PSI4 exiting successfully' in output
#            'MW' = { #OrderedDict of job molecular weights in g/mol
#                # job_id # : # molecular weight#
#                  '1'      : 12
#            }
#        }
#    }
#}
#    
#
#
#

##### NOTE: dft_methods is broken around line ~980 or so#####

##### ORIGINALLY ALL OF THESE WERE COMMENTED OUT. #####
##### I WILL BE UNCOMMENTING/REPLACING AS NEEDED. #####

import collections
import psi4
#import molutil
#import proc
import sys
import copy
#import functional
#import driver
# driver no longer exists, see proc.py and procrouting folder
#import atexit
#import p4const
import math
import os
#import p4util
#import re
# need json for new output structure
import json
# need numpy for holding property tensors
import numpy as np

def clean_up(db):
    try:
        # try to close the database
        db.close()
    except:
        # apparently it was already closed
        pass


## RESULTS LIST
#  A dictionary containing the quantities of interest for different levels of
#  theory and job_types
results = {
    'scf': {
        'scf_energy': 0,
        'scf_dipole': 0,
    },
    'cc2': {
        'correlation_energy': 0,
        'polarizability'    : 0,
        'optical_rotation'  : 0,
        'dipole_quadrupole_polarizability_tensor': 0,
    },
}

scalar_results = [
    'scf_energy',
    'cc_energy'
]

fscalar_results = [
    'polarizability',
    'rotation'
]

vector_results = [
    'scf_dipole'
]

tensor_results = [
    'quadrupole'
]

ftensor_results = [
    'polarizability_tensor',
    'rotation_tensor'
]


def initialize_database(database, kwargs):
    database['initialized']      = False
    database['inputs_generated'] = False
    database['jobs_complete']    = False
    database['results_computed'] = False
    database['distributed']      = True
    database['cutoff']           = False
    database['num_threads']      = False
    database['bsse']             = ''
    database['pcm']              = False
    if 'properties' in kwargs:
        database['n_body_func'] = psi4.properties
        database['omega'] = psi4.core.get_option('CCRESPONSE','OMEGA')
        database['omega'].pop()
    else:
        database['n_body_func'] = psi4.energy
    database['methods']          = {}


def extend_database(database, kwargs):
    for method,n_body_max in database['methods'].items():
        database[method] = collections.OrderedDict()
        database[method]['n_body_max'] = n_body_max
        # Assume we'll have scf energy and dipole
        # TODO: check dft energy printing and extend this
        database[method]['results'] = ['scf_energy','scf_dipole']
        # DFT methods
        if method == 'b3lyp':
            database[method]['results'].append('quadrupole')
        # Add correlation energy
        # Coupled Cluster methods
        if method in ['cc2','ccsd','cc3','eom-cc2','eom-ccsd']:
            database[method]['results'].append('cc_energy')
        # OTHER CORRELATED METHODS
        # Add properties
        if 'properties' in kwargs:
        # NOTE: DISABLED TENSOR HARVESTING FOR NOW
        # BUT ROTATIONS FOR G09 GET THEM AUTOMATICALLY
            if 'polarizability' in kwargs['properties']:
                database[method]['results'].append('polarizability')
#                database[method]['results'].append('polarizability_tensor')
            if 'rotation' in kwargs['properties']:
                database[method]['results'].append('rotation')
#                database[method]['results'].append('rotation_tensor')
            if 'quadrupole' in kwargs['properties']:
                database[method]['results'].append('quadrupole')
            # Need roa also!
        # Build list of job categories for iterating
        farm = []
        for n in range(1,n_body_max+1):
            # uncorrected MBE
            farm.append(n)
            # MBCP
            if database['bsse'] == 'mbcp' and n > 1:
                farm.append('{}-{}r'.format(n,1))
            # VMFC
            if database['bsse'] == 'vmfc':
                for m in range(1,n):
                    farm.append('{}-{}r'.format(n,m))
        database[method]['farm'] = farm
        for field in farm:
            database[method][field] = collections.OrderedDict()
            database[method][field]['num_jobs_complete'] = 0
            database[method][field]['total_num_jobs']    = 0
            database[method][field]['job_status'] = collections.OrderedDict()
            database[method][field]['MW'] = collections.OrderedDict()
            for result in database[method]['results']:
                database[method][field][result] = collections.OrderedDict()
                # @correction contains the sum of the cooked_data
                # i.e. the portion of the property due to only this n_body level
                database[method][field][result]['correction'] = 0
                # @result contains the total property at respective n_body level
                database[method][field][result][result] = 0
                database[method][field][result]['vmfc_correction'] = 0
                database[method][field][result]['vmfc_approximation'] = 0
                database[method][field][result]['mbcp_correction'] = 0
                database[method][field][result]['mbcp_approximation'] = 0
                database[method][field][result]['raw_data'] = collections.OrderedDict()
                database[method][field][result]['cooked_data'] = collections.OrderedDict()
        #for n in range(1,n_body_max+1):
        #    database[method][n] = collections.OrderedDict()
        #    database[method][n]['num_jobs_complete'] = 0
        #    database[method][n]['total_num_jobs']    = 0
        #    database[method][n]['job_status'] = collections.OrderedDict()
        #    for result in database[method]['results']:
        #        database[method][n][result] = collections.OrderedDict()
        #        # @correction contains the sum of the cooked_data
        #        # i.e. the portion of the property due to only this n_body level
        #        database[method][n][result]['correction'] = 0
        #        # @result contains the total property at respective n_body level
        #        database[method][n][result][result] = 0
        #        database[method][n][result]['raw_data'] = collections.OrderedDict()
        #        database[method][n][result]['cooked_data'] = collections.OrderedDict()
        #if database['bsse'] == 'mbcp':
        #    for n in range(2,n_body_max+1):
        #        ghost_jobs = '{}-{}r'.format(n,1)
        #        database[method][ghost_jobs] = collections.OrderedDict()
        #        database[method][ghost_jobs]['num_jobs_complete'] = 0
        #        database[method][ghost_jobs]['total_num_jobs']    = 0
        #        database[method][ghost_jobs]['job_status'] = collections.OrderedDict()
        #        for result in database[method]['results']:
        #            database[method][ghost_jobs][result] = collections.OrderedDict()
        #            # @correction contains the sum of the cooked_data
        #            # i.e. the portion of the property due to only this n_body level
        #            database[method][ghost_jobs][result]['correction'] = 0
        #            # @result contains the total property at respective n_body level
        #            database[method][ghost_jobs][result][result] = 0
        #            database[method][ghost_jobs][result]['raw_data'] = collections.OrderedDict()
        #            database[method][ghost_jobs][result]['cooked_data'] = collections.OrderedDict()
        #if database['bsse'] == 'vmfc':
        #    for n in range(2,n_body_max+1):
        #        for m in range(1, n):
        #            ghost_jobs = '{}-{}r'.format(n,m)
        #            database[method][ghost_jobs] = collections.OrderedDict()
        #            database[method][ghost_jobs]['num_jobs_complete'] = 0
        #            database[method][ghost_jobs]['total_num_jobs']    = 0
        #            database[method][ghost_jobs]['job_status'] = collections.OrderedDict()
        #            for result in database[method]['results']:
        #                database[method][ghost_jobs][result] = collections.OrderedDict()
        #                # @correction contains the sum of the cooked_data
        #                # i.e. the portion of the property due to only this n_body level
        #                database[method][ghost_jobs][result]['correction'] = 0
        #                # @result contains the total property at respective n_body level
        #                database[method][ghost_jobs][result][result] = 0
        #                database[method][ghost_jobs][result]['raw_data'] = collections.OrderedDict()
        #                database[method][ghost_jobs][result]['cooked_data'] = collections.OrderedDict()



def banner(db):
    psi4.print_stdout('\n')
    psi4.print_stdout('\t\t ------------------------------ \n')
    psi4.print_stdout('\t\t             N_Body             \n')
    psi4.print_stdout('\n')
    psi4.print_stdout('\t\t Input parameters \n')
    psi4.print_stdout('\t\t ---------------- \n')
    psi4.print_stdout('\t\t Job Type: {} \n'.format(db['n_body_func'].__name__))
    psi4.print_stdout('\t\t Methods:\n')
    for method in db['methods']:
        psi4.print_stdout('\t\t * {}:\t n <= '
                       '{}\n'.format(method,db['methods'][method]))
    psi4.print_stdout('\t\t Distributed: {} \n'.format(db['distributed']))
    psi4.print_stdout('\n')
    # STATUS
    psi4.print_stdout('\t\t Current Status \n')
    psi4.print_stdout('\t\t -------------- \n')
    input_message = """\n
    Input files have been generated for the above listed methods. You
    should run these in whatever way you prefer. Invoking psi4 while they run
    will check the individual jobs status. As the jobs required for each level
    of approximation complete the data for the results listed above will be
    automatically collected and associated corrections computed. These will
    appear in the results section below as they become available.

    !! There are ancilliary python scripts that interact with the n_body
    database file for convenient submission of the generated input files and
    custom result manipulations that you may wish to perform. They are
    available here: github.com/bgpeyton/n_body\n"""
    complete_message = """\n
    All requested jobs have completed running. Results and corrections should be
    printed below for all approximation levels requested. If not try running
    psi4 one more time. If after analyzing the data you desire another method or
    higher approximation levels modify the input file appropriately and run
    psi4. The existing jobs will remain and input files for the newly requested
    jobs will appear. Upon submission and completion of these jobs the results
    section below will show the updated results in addition to those currently
    shown.\n"""
    # Input files have been generated
    if db['inputs_generated']:
        psi4.print_stdout(input_message)
    # All the jobs have completed
    if db['jobs_complete']:
        psi4.print_stdout(complete_message)
    else:
        for method in db['methods']:
            psi4.print_stdout('        ==> {} jobs running <==\n'.format(method.upper()))
            # more super annoying prints, muting for now
#            for n in range(1,db[method]['n_body_max']+1):
#                for job,stat in db[method][n]['job_status'].items():
#                    if stat != 'complete':
#                        psi4.print_stdout('{} {}\n'.format(job,stat))

#    # RESULTS
#    #### NOTE: SKIPPING RESULTS FOR NOW #####
#    psi4.print_stdout('\t\t Results \n')
#    psi4.print_stdout('\t\t ------- \n')
#    for method in db['methods']:
#        psi4.print_stdout('        ==> {} <==\n'.format(method.upper()))
#        if psi4.get_global_option('PRINT') < 2:
#            psi4.print_stdout('    Tensor quantities not being printed. To print them '
#                           'increase print > 1 and rerun psi4.\n')
#        for result in db[method]['results']:
#            for n in range(1,db[method]['n_body_max']+1):
#                # Does this quantity exist?
#                if db[method][n][result][result]:
#                    # What type of quantity is it??
#                    if result in scalar_results:
#                        print_scalar(db, method, n, result)
#                    elif result in vector_results:
#                        print_vector(db, method, n, result)
#                    elif result in fscalar_results:
#                        print_fscalar(db, method, n, result)
#                    elif len(db[method][n][result][result]) % 9 == 0:
#                        print_rank_two(db, method, n, result)
#                    else:
#                        getattr(sys.modules[__name__],'print_{}_results'.format(result))(db,method,n)
#        #if 'rotation' in db[method]['results']:
#        #    for n in range(1,db[method]['n_body_max']+1):
#        #        print_mass_adjusted_rotation(db, method, n)
    psi4.print_stdout('\t\t ------------------------------ \n')

def print_scalar(db, method, n, result):
    # Print header
    if n == 1:
        psi4.print_stdout('{:<58}'.format('        {}:    '.format(re.sub('_',' ',result).upper())))
        psi4.print_stdout('\n\n')
    # Print scalars
    corr = db[method][n][result]['correction']
    approx = db[method][n][result][result]
    psi4.print_stdout('{:<58}'.format('            {} Correction:'.format(print_body(n))))
    psi4.print_stdout('    {:> 22.14f}\n'.format(corr[0]))
    psi4.print_stdout('{:<58}'.format('            {} Approximation:'.format(print_body(n))))
    psi4.print_stdout('    {:> 22.14f}\n'.format(approx[0]))
    if db['bsse'] == ('vmfc' or 'mbcp') and n > 1:
        mbcp_corr = db[method][n][result]['mbcp_correction']
        mbcp_approx = db[method][n][result]['mbcp_approximation']
        psi4.print_stdout('{:<58}'.format('            MBCP({}) Correction:'.format(n)))
        psi4.print_stdout('    {:> 22.14f}\n'.format(mbcp_corr[0]))
        psi4.print_stdout('{:<58}'.format('            MBCP({}) Approximation:'.format(n)))
        psi4.print_stdout('    {:> 22.14f}\n'.format(mbcp_approx[0]))
    if db['bsse'] == 'vmfc' and n > 1:
        vmfc_corr = db[method][n][result]['vmfc_correction']
        vmfc_approx = db[method][n][result]['vmfc_approximation']
        psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Correction:'.format(print_body(n))))
        psi4.print_stdout('    {:> 22.14f}\n'.format(vmfc_corr[0]))
        psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Approximation:'.format(print_body(n))))
        psi4.print_stdout('    {:> 22.14f}\n'.format(vmfc_approx[0]))
        psi4.print_stdout('{:<58}'.format('            VM {0} - {0} Approximation:'.format(print_body(n))))
        psi4.print_stdout('    {:> 22.14f}\n'.format(vmfc_approx[0] - approx[0]))
    # Print errors (but only if all fragments have been computed)
    molecule = psi4.get_active_molecule()
    if db[method].has_key(molecule.nfragments()):
        if db[method][molecule.nfragments()][result][result]:
            actual = db[method][molecule.nfragments()][result][result]
            error = [x-y for x,y in zip(approx,actual)]
            p_error = [(x/abs(y))*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            {} Error:'.format(print_body(n))))
            psi4.print_stdout('    {:> 22.14f}\n'.format(error[0]))
            psi4.print_stdout('{:<58}'.format('            {} Percent Error:'.format(print_body(n))))
            psi4.print_stdout('    {:> 22.14f}\n'.format(p_error[0]))
        if db['bsse'] == ('vmfc' or 'mbcp') and n > 1:
            mbcp_corr = copy.deepcopy(db[method][n][result]['mbcp_correction'])
            mbcp_approx = copy.deepcopy(db[method][n][result]['mbcp_approximation'])
            actual = db[method][molecule.nfragments()][result]['mbcp_approximation']
            error = [x-y for x,y in zip(mbcp_approx,actual)]
            p_error = [(x/abs(y))*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            MBCP({}) Error:'.format(n)))
            while error:
                psi4.print_stdout('    {:> 22.14f}'.format(error.pop(0)))
            psi4.print_stdout('\n')
            psi4.print_stdout('{:<58}'.format('            MBCP({}) Percent Error:'.format(n)))
            while p_error:
                psi4.print_stdout('    {:> 22.14f}'.format(p_error.pop(0)))
            psi4.print_stdout('\n')
        if db['bsse'] == 'vmfc' and n > 1:
            vmfc_corr = copy.deepcopy(db[method][n][result]['vmfc_correction'])
            vmfc_approx = copy.deepcopy(db[method][n][result]['vmfc_approximation'])
            actual = db[method][molecule.nfragments()][result]['vmfc_approximation']
            error = [x-y for x,y in zip(vmfc_approx,actual)]
            p_error = [(x/abs(y))*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Error:'.format(print_body(n))))
            while error:
                psi4.print_stdout('    {:> 22.14f}'.format(error.pop(0)))
            psi4.print_stdout('\n')
            psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Percent Error:'.format(print_body(n))))
            while p_error:
                psi4.print_stdout('    {:> 22.14f}'.format(p_error.pop(0)))
            psi4.print_stdout('\n')
    # Print Interaction scalars 
    if n > 1:
        mono_approx = db[method][1][result][result]
        interact = [x-y for x,y in zip(approx,mono_approx)]
        psi4.print_stdout('{:<58}'.format('            {} Interaction Energy Approximation:'.format(print_body(n))))
        psi4.print_stdout('    {:> 22.14f}\n'.format(interact[0]))
        # Print the Interaction errors
        if db[method].has_key(molecule.nfragments()):
            if db[method][molecule.nfragments()][result][result]:
                actual = db[method][molecule.nfragments()][result][result]
                actual_interact = [x-y for x,y in zip(actual,mono_approx)]
                error = [x-y for x,y in zip(interact,actual_interact)]
                p_error = [(x/abs(y))*100 for x,y in zip(error,actual_interact)]
                psi4.print_stdout('{:<58}'.format('            {} Interaction Error:'.format(print_body(n))))
                psi4.print_stdout('    {:> 22.14f}\n'.format(error[0]))
                psi4.print_stdout('{:<58}'.format('            {} Interaction Percent Error:'.format(print_body(n))))
                psi4.print_stdout('    {:> 22.14f}\n\n'.format(p_error[0]))
        else:
            psi4.print_stdout('\n')
    else:
        psi4.print_stdout('\n')

def print_vector(db, method, n, result):
    # Print header
    if n == 1:
        print(result)
        psi4.print_stdout('{:<58}'.format('        {}:    '.format(re.sub('_',' ',result).upper())))
        psi4.print_stdout('    {:>14}'.format('X-component'))
        psi4.print_stdout('    {:>14}'.format('Y-component'))
        psi4.print_stdout('    {:>14}'.format('Z-component'))
        psi4.print_stdout('    {:>14}'.format('Total'))
        psi4.print_stdout('\n\n')
    # Print components
    if db[method][n][result][result]:
        corr = copy.deepcopy(db[method][n][result]['correction']) 
        approx = copy.deepcopy(db[method][n][result][result])
        psi4.print_stdout('{:<58}'.format('            {} Correction:'.format(print_body(n))))
        total = 0
        while corr:
            current = corr.pop(0)
            psi4.print_stdout('    {:> 14.6f}'.format(current))
            total += current ** 2
        psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(total)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            {} Approximation:'.format(print_body(n))))
        total = 0
        while approx:
            current = approx.pop(0)
            psi4.print_stdout('    {:> 14.6f}'.format(current))
            total += current ** 2
        psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(total)))
        psi4.print_stdout('\n')
    if db['bsse'] == ('vmfc' or 'mbcp') and n > 1:
        mbcp_corr = copy.deepcopy(db[method][n][result]['mbcp_correction'])
        mbcp_approx = copy.deepcopy(db[method][n][result]['mbcp_approximation'])
        psi4.print_stdout('{:<58}'.format('            MBCP({}) Correction:'.format(n)))
        total = 0
        while mbcp_corr:
            current = mbcp_corr.pop(0)
            psi4.print_stdout('    {:> 14.6f}'.format(current))
            total += current ** 2
        psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(total)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            MBCP({}) Approximation:'.format(n)))
        total = 0
        while mbcp_approx:
            current = mbcp_approx.pop(0)
            psi4.print_stdout('    {:> 14.6f}'.format(current))
            total += current ** 2
        psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(total)))
        psi4.print_stdout('\n')
    if db['bsse'] == 'vmfc' and n > 1:
        vmfc_corr = copy.deepcopy(db[method][n][result]['vmfc_correction'])
        vmfc_approx = copy.deepcopy(db[method][n][result]['vmfc_approximation'])
        psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Correction:'.format(print_body(n))))
        total = 0
        while vmfc_corr:
            current = vmfc_corr.pop(0)
            psi4.print_stdout('    {:> 14.6f}'.format(current))
            total += current ** 2
        psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(total)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Approximation:'.format(print_body(n))))
        total = 0
        while vmfc_approx:
            current = vmfc_approx.pop(0)
            psi4.print_stdout('    {:> 14.6f}'.format(current))
            total += current ** 2
        psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(total)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            VM {0} - {0} Approximation:'.format(print_body(n))))
        approx = copy.deepcopy(db[method][n][result][result])
        vmfc_approx = copy.deepcopy(db[method][n][result]['vmfc_approximation'])
        total = 0
        while vmfc_approx:
            current = approx.pop(0)
            vmfc_current = vmfc_approx.pop(0)
            psi4.print_stdout('    {:> 14.6f}'.format(vmfc_current - current))
            total += (vmfc_current - current) ** 2
        psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(total)))
        psi4.print_stdout('\n')
    # Print errors (but only if all fragments have been computed)
    #
    # ERROR printing is accomplished by computing an error vector. This vector takes into account both the
    # direction and magnitude errors. The total error is reported as the magnitude of this error vector and
    # not the error in the magnitude of the approximate vector. This reported error will always be larger 
    # than the error in magnitude between the approximate and true vectors as it also incorporates directional
    # error.
    #
    # PERCENT ERROR printing is slightly different as creating a percent error vector makes no sense. Instead
    # the percent error of each component is calculated separately, as the error component divided by the true component
    # including sign (which denotes in which direction the approximate vector is from the true vector) times a
    # hundred. But the total percent error is reported as the magnitude of the error vector divided by the magnitude of
    # the true vector, not the magnitude of the "percent error vector".
    molecule = psi4.get_active_molecule()
    if db[method].has_key(molecule.nfragments()):
        if db[method][molecule.nfragments()][result][result]:
            approx = copy.deepcopy(db[method][n][result][result])
            actual = db[method][molecule.nfragments()][result][result]
            actual_total = math.sqrt(sum([x*y for x,y in zip(actual,actual)]))
            error = [x-y for x,y in zip(approx,actual)]
            p_error = [(x/y)*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            {} Error:'.format(print_body(n))))
            error_total = 0
            while error:
                current = error.pop(0)
                psi4.print_stdout('    {:> 14.6f}'.format(current))
                error_total += current ** 2
            psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(error_total)))
            psi4.print_stdout('\n')
            psi4.print_stdout('{:<58}'.format('            {} Percent Error:'.format(print_body(n))))
            while p_error:
                current = p_error.pop(0)
                psi4.print_stdout('    {:> 14.6f}'.format(current))
            psi4.print_stdout('    {:> 14.6f}'.format((math.sqrt(error_total) / actual_total) * 100))
            psi4.print_stdout('\n')
        if db['bsse'] == ('vmfc' or 'mbcp') and n > 1:
            mbcp_corr = copy.deepcopy(db[method][n][result]['mbcp_correction'])
            mbcp_approx = copy.deepcopy(db[method][n][result]['mbcp_approximation'])
            actual = db[method][molecule.nfragments()][result]['mbcp_approximation']
            actual_total = math.sqrt(sum([x*y for x,y in zip(actual,actual)]))
            error = [x-y for x,y in zip(mbcp_approx,actual)]
            p_error = [(x/abs(y))*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            MBCP({}) Error:'.format(n)))
            error_total = 0
            while error:
                current = error.pop(0)
                psi4.print_stdout('    {:> 14.6f}'.format(current))
                error_total += current ** 2
            psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(error_total)))
            psi4.print_stdout('\n')
            psi4.print_stdout('{:<58}'.format('            MBCP({}) Percent Error:'.format(n)))
            while p_error:
                current = p_error.pop(0)
                psi4.print_stdout('    {:> 14.6f}'.format(current))
            psi4.print_stdout('    {:> 14.6f}'.format((math.sqrt(error_total) / actual_total) * 100))
            psi4.print_stdout('\n')
        if db['bsse'] == 'vmfc' and n > 1:
            vmfc_corr = copy.deepcopy(db[method][n][result]['vmfc_correction'])
            vmfc_approx = copy.deepcopy(db[method][n][result]['vmfc_approximation'])
            actual = db[method][molecule.nfragments()][result]['vmfc_approximation']
            actual_total = math.sqrt(sum([x*y for x,y in zip(actual,actual)]))
            error = [x-y for x,y in zip(vmfc_approx,actual)]
            p_error = [(x/abs(y))*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Error:'.format(print_body(n))))
            while error:
                current = error.pop(0)
                psi4.print_stdout('    {:> 14.6f}'.format(current))
                error_total += current ** 2
            psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(error_total)))
            psi4.print_stdout('\n')
            psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Percent Error:'.format(print_body(n))))
            while p_error:
                current = p_error.pop(0)
                psi4.print_stdout('    {:> 14.6f}'.format(current))
            psi4.print_stdout('    {:> 14.6f}'.format((math.sqrt(error_total) / actual_total) * 100))
            psi4.print_stdout('\n')
    # Print Interaction properties (Is this even meaningful?)
    if n > 1:
        approx = copy.deepcopy(db[method][n][result][result])
        mono_approx = db[method][1][result][result]
        interact = [x-y for x,y in zip(approx,mono_approx)]
        psi4.print_stdout('{:<58}'.format('            {} Interaction Approximation:'.format(print_body(n))))
        total = 0
        while interact: 
            current = interact.pop(0)
            psi4.print_stdout('    {:> 14.6f}'.format(current))
            total += current ** 2
        psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(total)))
        psi4.print_stdout('\n')
        # Print the Interaction energy errors
        if db[method].has_key(molecule.nfragments()):
            if db[method][molecule.nfragments()][result][result]:
                actual = db[method][molecule.nfragments()][result][result]
                actual_interact = [x-y for x,y in zip(actual,mono_approx)]
                actual_inter_total = math.sqrt(sum([x*y for x,y in zip(actual_interact,actual_interact)]))
                interact = [x-y for x,y in zip(approx,mono_approx)]
                error = [x-y for x,y in zip(interact,actual_interact)]
                p_error = [(x/y)*100 for x,y in zip(error,actual_interact)]
                psi4.print_stdout('{:<58}'.format('            {} Interaction Error:'.format(print_body(n))))
                error_total = 0
                while error:
                    current = error.pop(0)
                    psi4.print_stdout('    {:> 14.6f}'.format(current))
                    error_total += current ** 2
                psi4.print_stdout('    {:> 14.6f}'.format(math.sqrt(error_total)))
                psi4.print_stdout('\n')
                psi4.print_stdout('{:<58}'.format('            {} Interaction Percent Error:'.format(print_body(n))))
                while p_error:
                    current = p_error.pop(0)
                    psi4.print_stdout('    {:> 14.6f}'.format(current))
                psi4.print_stdout('    {:> 14.6f}'.format((math.sqrt(error_total) / actual_inter_total) * 100))
                psi4.print_stdout('\n\n')
        else:
            psi4.print_stdout('\n')
    else:
        psi4.print_stdout('\n')

def print_fscalar(db, method, n, result):
#def print_rotation_results(db, method, n):
    omega = psi4.get_global_option('OMEGA')
    if len(omega) > 1:
        units = omega.pop()
    # Print header
    if n == 1:
        psi4.print_stdout('{:<58}'.format('        {}:    '.format(re.sub('_',' ',result).upper())))
        while omega:
            psi4.print_stdout('    {:> 9.2f} {}'.format(omega.pop(0),units))
        psi4.print_stdout('\n\n')
    # Print n-body results
    if db[method][n][result][result]:
        corr = copy.deepcopy(db[method][n][result]['correction']) 
        approx = copy.deepcopy(db[method][n][result][result])
        psi4.print_stdout('{:<58}'.format('            {} Correction:'.format(print_body(n))))
        while corr:
            psi4.print_stdout('    {:> 12.6f}'.format(corr.pop(0)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            {} Approximation:'.format(print_body(n))))
        while approx:
            psi4.print_stdout('    {:> 12.6f}'.format(approx.pop(0)))
        #psi4.print_stdout('\n')
        psi4.print_stdout('\n            {:->118}\n'.format(''))
        psi4.print_stdout('{:<58}\n'.format('            Cluster Contributions: (Distance)'))
        # Print cluster contributions
        mol = psi4.get_active_molecule()
        indexes = molutil.extract_cluster_indexing(mol, n)
        for index in indexes:
            index_dir = proc.cluster_dir(index)
            distance  = cluster_distance(n,index)
            try:
                psi4.print_stdout('{:<58}'.format('              {}                  ({: 8.2f})'.format(index_dir,distance)))
            except:
                psi4.print_stdout('{:<58}'.format('              {}                           '.format(index_dir)))
            current = db[method][n][result]['cooked_data'][index_dir]
            while current:
                psi4.print_stdout('    {:> 12.6f}'.format(current.pop(0)))
            psi4.print_stdout('\n')
        psi4.print_stdout('\n            {:->118}\n'.format(''))

    # Print BSSE results
    if db['bsse'] == ('vmfc' or 'mbcp') and n > 1:
        mbcp_corr = copy.deepcopy(db[method][n][result]['mbcp_correction'])
        mbcp_approx = copy.deepcopy(db[method][n][result]['mbcp_approximation'])
        psi4.print_stdout('{:<58}'.format('            MBCP({}) Correction:'.format(n)))
        while mbcp_corr:
            psi4.print_stdout('    {:> 12.6f}'.format(mbcp_corr.pop(0)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            MBCP({}) Approximation:'.format(n)))
        while mbcp_approx:
            psi4.print_stdout('    {:> 12.6f}'.format(mbcp_approx.pop(0)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            MBCP({0}) - {1} Approximation:'.format(n,print_body(n))))
        approx = copy.deepcopy(db[method][n][result][result])
        mbcp_approx = copy.deepcopy(db[method][n][result]['mbcp_approximation'])
        while mbcp_approx:
            psi4.print_stdout('    {:> 12.6f}'.format(mbcp_approx.pop(0) - approx.pop(0)))
        psi4.print_stdout('\n')
    if db['bsse'] == 'vmfc' and n > 1:
        vmfc_corr = copy.deepcopy(db[method][n][result]['vmfc_correction'])
        vmfc_approx = copy.deepcopy(db[method][n][result]['vmfc_approximation'])
        psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Correction:'.format(print_body(n))))
        while vmfc_corr:
            psi4.print_stdout('    {:> 12.6f}'.format(vmfc_corr.pop(0)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Approximation:'.format(print_body(n))))
        while vmfc_approx:
            psi4.print_stdout('    {:> 12.6f}'.format(vmfc_approx.pop(0)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            VM {0} - {0} Approximation:'.format(print_body(n))))
        approx = copy.deepcopy(db[method][n][result][result])
        vmfc_approx = copy.deepcopy(db[method][n][result]['vmfc_approximation'])
        while vmfc_approx:
            psi4.print_stdout('    {:> 12.6f}'.format(vmfc_approx.pop(0) - approx.pop(0)))
        psi4.print_stdout('\n')
    # Print errors (but only if all fragments have been computed)
    molecule = psi4.get_active_molecule()
    if db[method].has_key(molecule.nfragments()):
        if db[method][molecule.nfragments()][result][result]:
            approx = copy.deepcopy(db[method][n][result][result])
            actual = db[method][molecule.nfragments()][result][result]
            error = [x-y for x,y in zip(approx,actual)]
            p_error = [(x/abs(y))*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            {} Error:'.format(print_body(n))))
            while error:
                psi4.print_stdout('    {:> 12.6f}'.format(error.pop(0)))
            psi4.print_stdout('\n')
            psi4.print_stdout('{:<58}'.format('            {} Percent Error:'.format(print_body(n))))
            while p_error:
                psi4.print_stdout('    {:> 12.6f}'.format(p_error.pop(0)))
            psi4.print_stdout('\n')
        if db['bsse'] == ('vmfc' or 'mbcp') and n > 1:
            mbcp_corr = copy.deepcopy(db[method][n][result]['mbcp_correction'])
            mbcp_approx = copy.deepcopy(db[method][n][result]['mbcp_approximation'])
            actual = db[method][molecule.nfragments()][result]['mbcp_approximation']
            error = [x-y for x,y in zip(mbcp_approx,actual)]
            p_error = [(x/abs(y))*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            MBCP({}) Error:'.format(n)))
            while error:
                psi4.print_stdout('    {:> 12.6f}'.format(error.pop(0)))
            psi4.print_stdout('\n')
            psi4.print_stdout('{:<58}'.format('            MBCP({}) Percent Error:'.format(n)))
            while p_error:
                psi4.print_stdout('    {:> 12.6f}'.format(p_error.pop(0)))
            psi4.print_stdout('\n')
        if db['bsse'] == 'vmfc' and n > 1:
            vmfc_corr = copy.deepcopy(db[method][n][result]['vmfc_correction'])
            vmfc_approx = copy.deepcopy(db[method][n][result]['vmfc_approximation'])
            actual = db[method][molecule.nfragments()][result]['vmfc_approximation']
            error = [x-y for x,y in zip(vmfc_approx,actual)]
            p_error = [(x/abs(y))*100 for x,y in zip(error,actual)]
            psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Error:'.format(print_body(n))))
            while error:
                psi4.print_stdout('    {:> 12.6f}'.format(error.pop(0)))
            psi4.print_stdout('\n')
            psi4.print_stdout('{:<58}'.format('            Valiron-Mayer {} Percent Error:'.format(print_body(n))))
            while p_error:
                psi4.print_stdout('    {:> 12.6f}'.format(p_error.pop(0)))
            psi4.print_stdout('\n')
    # Print Interaction properties (Is this even meaningful?)
    if n > 1:
        mono_approx = db[method][1][result][result]
        interact = [x-y for x,y in zip(approx,mono_approx)]
        psi4.print_stdout('{:<58}'.format('            {} Interaction Approximation:'.format(print_body(n))))
        while interact: 
            psi4.print_stdout('    {:> 12.6f}'.format(interact.pop(0)))
        psi4.print_stdout('\n')
        # Print the Interaction energy errors
        if db[method].has_key(molecule.nfragments()):
            if db[method][molecule.nfragments()][result][result]:
                actual = db[method][molecule.nfragments()][result][result]
                actual_interact = [x-y for x,y in zip(actual,mono_approx)]
                interact = [x-y for x,y in zip(approx,mono_approx)]
                error = [x-y for x,y in zip(interact,actual_interact)]
                p_error = [(x/abs(y))*100 for x,y in zip(error,actual_interact)]
                psi4.print_stdout('{:<58}'.format('            {} Interaction Error:'.format(print_body(n))))
                while error:
                    psi4.print_stdout('    {:> 12.6f}'.format(error.pop(0)))
                psi4.print_stdout('\n')
                psi4.print_stdout('{:<58}'.format('            {} Interaction Percent Error:'.format(print_body(n))))
                while p_error:
                    psi4.print_stdout('    {:> 12.6f}'.format(p_error.pop(0)))
                psi4.print_stdout('\n\n')
        else:
            psi4.print_stdout('\n')
    else:
        psi4.print_stdout('\n')

def print_cc_energy_results(db, method, n):
    if n == 1:
        psi4.print_stdout('{:<58}'.format('        CC Energy:     '))
        psi4.print_stdout('\n\n')
    # Print CC energies
    if db[method][n]['cc_energy']['cc_energy']:
        corr = db[method][n]['cc_energy']['correction']
        approx = db[method][n]['cc_energy']['cc_energy']
        psi4.print_stdout('{:<58}'.format('            {} Correction:'.format(print_body(n))))
        psi4.print_stdout('    {:> 22.14f}\n'.format(corr[0]))
        psi4.print_stdout('{:<58}'.format('            {} Approximation:'.format(print_body(n))))
        psi4.print_stdout('    {:> 22.14f}\n'.format(approx[0]))
        # Print CC energy errors
        #molecule = psi4.get_active_molecule()
        #if db[method][molecule.nfragments()]['cc_energy']['']
        # Print CC Interaction energies
        # Print CC Interaction energies errors
        if n > 1:
            mono_approx = db[method][1]['cc_energy']['cc_energy']
            interact = [x-y for x,y in zip(approx,mono_approx)]
            psi4.print_stdout('{:<58}'.format('            {} Interaction Energy Approximation:'.format(print_body(n))))
            psi4.print_stdout('    {:> 22.14f}\n\n'.format(interact[0]))
        else:
            psi4.print_stdout('\n')


def print_polarizability_results(db, method, n):
    omega = psi4.get_global_option('OMEGA')
    if len(omega) > 1:
        units = omega.pop()
    if n == 1:
        psi4.print_stdout('{:<58}'.format('        POLARIZABILITY:    '))
        while omega:
            psi4.print_stdout('    {:> 9.2f} {}'.format(omega.pop(0),units))
        psi4.print_stdout('\n\n')
    if db[method][n]['polarizability']['polarizability']:
        corr = copy.deepcopy(db[method][n]['polarizability']['correction']) 
        approx = copy.deepcopy(db[method][n]['polarizability']['polarizability'])
        psi4.print_stdout('{:<58}'.format('            {} Correction:'.format(print_body(n))))
        while corr:
            psi4.print_stdout('    {:> 12.6f}'.format(corr.pop(0)))
        psi4.print_stdout('\n')
        psi4.print_stdout('{:<58}'.format('            {} Approximation:'.format(print_body(n))))
        while approx:
            psi4.print_stdout('    {:> 12.6f}'.format(approx.pop(0)))
        psi4.print_stdout('\n\n')

def print_polarizability_tensor_results(db, method, n):
    if db[method][n]['polarizability_tensor']['polarizability_tensor']:
        corr = copy.deepcopy(db[method][n]['polarizability_tensor']['correction'])
        approx = copy.deepcopy(db[method][n]['polarizability_tensor']['polarizability_tensor'])
        if psi4.get_global_option('PRINT') > 1:
            psi4.print_stdout('        {} Correction:   \n'.format(print_body(n)))
            while corr:
                psi4.print_stdout('       {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n'.format(corr.pop(0)))
                psi4.print_stdout('       {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n'.format(corr.pop(0)))
                psi4.print_stdout('       {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n\n'.format(corr.pop(0)))

            psi4.print_stdout('        {} Approximation: \n'.format(print_body(n)))
            while approx:
                psi4.print_stdout('       {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n'.format(approx.pop(0)))
                psi4.print_stdout('       {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n'.format(approx.pop(0)))
                psi4.print_stdout('       {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n\n'.format(approx.pop(0)))

def print_rank_two(db, method, n, result):
    print('Rank two quantity: {}'.format(result))
 

def mass_adjust_rotations(db, method, n):
    #omega = psi4.get_global_option('OMEGA')
    #if len(omega) > 1:
    #    units = omega.pop()
    # Print result header
    #if n == 1:
    #    psi4.print_stdout('{:<58}'.format('        MASS-ADJUSTED OPTICAL ROTATION:    '))
    #    while omega:
    #        psi4.print_stdout('    {:> 11.2f} {}'.format(omega.pop(0),units))
    #    psi4.print_stdout('\n\n')

    # Compute conversion factors
    psi_pi = 3.14159265358979323846264338327950288
    bohr2a4 = p4const.psi_bohr2angstroms ** 4
    m2a = p4const.psi_bohr2angstroms * 1.0e-10
    hbar = p4const.psi_h/(2.0 * psi_pi);
    prefactor = 1.0e-2 * hbar/(p4const.psi_c * 2.0 * psi_pi * p4const.psi_me * m2a * m2a);
    prefactor *= prefactor;
    prefactor *= 288.0e-30 * psi_pi * psi_pi * p4const.psi_na * bohr2a4;

    # Get omega list and convert to atomic units
    omega = psi4.get_global_option('OMEGA')
    units = 'AU'
    if len(omega) > 1:
        units = omega.pop()
    if units in ['HZ', 'Hz', 'hz']:
        omega = [x * (p4const.psi_h / p4const.psi_hartree2J) for x in omega]
    elif units in ['AU', 'Au', 'au']:
        pass
    elif units in ['NM', 'nm']:
        omega = [((p4const.psi_c * p4const.psi_h * 1e9) / (x * p4const.psi_hartree2J)) for x in omega]
    elif units in ['EV', 'ev', 'eV']:
        omega = [(x / p4const.psi_hartree2ev) for x in omega]

    # Molecular mass (assuming frag 1 is solute)
    mol = psi4.get_active_molecule()
    solute = mol.extract_subsets(1)
    M = 0
    for k in range(solute.natom()):
        M += solute.mass(k)

    # Conversion factor (tensor to specific rotation)
    ### G09 ###
    g09_conv = -2155513.07465937873759239624
    ### PSI4 ###
    # Conversion factor only works for MVG data
    psi4_conv = prefactor
    psi4_conv /= 3.0

    if method in ['cc2','ccsd','cc3','eom-cc2','eom-ccsd']:
        conv = psi4_conv
    else:
        conv = g09_conv

    # Print results
    if db[method][n]['rotation_tensor']['rotation_tensor']:
        # Make a list of what to mass adjust and where to store them
        mass_adjustees = ['correction', 'rotation_tensor']
        storees = ['correction', 'rotation']
        if db['bsse'] == 'vmfc' and n > 1:
            mass_adjustees.extend(['vmfc_correction', 'vmfc_approximation', 'mbcp_correction', 'mbcp_approximation'])
            storees.extend(['vmfc_correction', 'vmfc_approximation', 'mbcp_correction', 'mbcp_approximation'])
        elif db['bsse'] == 'mbcp' and n > 1:
            mass_adjustees.extend(['mbcp_correction', 'mbcp_approximation'])
            storees.extend(['mbcp_correction', 'mbcp_approximation'])

        # Compute mass adjusted rotation corrections
        for adjustee in mass_adjustees:
            mass_adjusted_rotations = []
            omega_copy = copy.deepcopy(omega)
            data = copy.deepcopy(db[method][n]['rotation_tensor'][adjustee])
            while data:
                trace = 0
                # Tensor trace
                trace += data.pop(0)  # Element 1
                data.pop(0)           # Element 2
                data.pop(0)           # Element 3
                data.pop(0)           # Element 4
                trace += data.pop(0)  # Element 5
                data.pop(0)           # Element 6
                data.pop(0)           # Element 7
                data.pop(0)           # Element 8
                trace += data.pop(0)  # Element 9
                nu = omega_copy.pop(0)
                # Convert to specific rotation
                #ctrace *= -2155513.07465937873759239624
                trace *= conv
                if method not in ['cc2','ccsd','cc3','eom-cc2','eom-ccsd']:
                    trace *= nu * nu
                trace /= M
                mass_adjusted_rotations.append(trace)
            # Overwrite rotation corrections with mass adjusted rotation corrections 
            #storee = mass_adjusted_rotations
            storee = storees.pop(0)
            db[method][n]['rotation'][storee] = mass_adjusted_rotations

        # Compute individual pair contributions
        # Future if block that turns on/off pair-analysis printing
        ##psi4.print_stdout('\n            {:->118}\n'.format(''))
        ##psi4.print_stdout('{:<58}\n'.format('            Cluster Contributions: (Distance)'))
        indexes = molutil.extract_cluster_indexing(mol, n)
        for index in indexes:
            imass_adjusted_rotations = []
            index_dir = proc.cluster_dir(index)
            #distance = cluster_distance(n,index)
            #try:
            #    psi4.print_stdout('{:<58}'.format('              {}                 ({: 6.2f})'.format(index_dir,distance)))
            #except:
            #    psi4.print_stdout('{:<58}'.format('              {}                           '.format(index_dir)))
            individual = copy.deepcopy(db[method][n]['rotation_tensor']['cooked_data'][index_dir])
            omega_copy = copy.deepcopy(omega)
            while individual:
                itrace = 0
                # Tensor trace
                itrace += individual.pop(0) # Element 1
                individual.pop(0)           # Element 2
                individual.pop(0)           # Element 3
                individual.pop(0)           # Element 4
                itrace += individual.pop(0) # Element 5
                individual.pop(0)           # Element 6
                individual.pop(0)           # Element 7
                individual.pop(0)           # Element 8
                itrace += individual.pop(0) # Element 9
                nu = omega_copy.pop(0)
                # Convert to specific rotation
                #itrace *= -2155513.07465937873759239624
                itrace *= conv
                if method not in ['cc2','ccsd','cc3','eom-cc2','eom-ccsd']:
                    itrace *= nu * nu
                itrace /= M
                imass_adjusted_rotations.append(itrace)
            # Overwrite rotation with mass-adjusted rotation
            db[method][n]['rotation']['cooked_data'][index_dir] = imass_adjusted_rotations

def print_rotation_tensor_results(db, method, n):
    if db[method][n]['rotation_tensor']['rotation_tensor']:
        corr = copy.deepcopy(db[method][n]['rotation_tensor']['correction'])
        approx = copy.deepcopy(db[method][n]['rotation_tensor']['rotation_tensor'])
        if psi4.get_global_option('PRINT') > 1:
            psi4.print_stdout('        {} Correction:   \n'.format(print_body(n)))
            while corr:
                psi4.print_stdout('       {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n'.format(corr.pop(0)))
                psi4.print_stdout('       {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n'.format(corr.pop(0)))
                psi4.print_stdout('       {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(corr.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n\n'.format(corr.pop(0)))

            psi4.print_stdout('        {} Approximation: \n'.format(print_body(n)))
            while approx:
                psi4.print_stdout('       {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n'.format(approx.pop(0)))
                psi4.print_stdout('       {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n'.format(approx.pop(0)))
                psi4.print_stdout('       {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}'.format(approx.pop(0)))
                psi4.print_stdout('    {: 10.6f}\n\n'.format(approx.pop(0)))

def print_body(n):
    num_2_word = ['Zero', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven',
                  'Eight', 'Nine', 'Ten', 'Eleven', 'Twelve', 'Thirteen']

    return('{}-Body'.format(num_2_word[n]))

    

##### NOTE: BROKEN #####
##### AJ SAYS IT NO LONGER EXISTS #####
##### TALK TO HIM #####
# Build list of dft functional names
#dft_methods = []
#for ssuper in functional.superfunctional_list():
#    dft_methods.append(ssuper.name().lower())
## Adds HF for now as psi4 can't yet do optical rotations
#dft_methods.append('hf')


def process_options(name, db, options):
    processed_options = {'methods':{}}
    molecule = psi4.core.get_active_molecule()
    nfrags = molecule.nfragments()

    try:
        if isinstance(options, bool):
        # Assume user desires to run jobs as 'name' and n_body_max = nfrags
            if options: # n_body = True in input
                processed_options['methods'].update({name.lower():nfrags})
            else:
                raise Exception('n_body = False is invalid. If you meant to run '
                                'without n_body, just remove it from the input.')
        elif isinstance(options, int):
        # Assume user desires to run jobs as 'name' and n_body_max = options
            processed_options['methods'].update({name.lower():options})
        elif isinstance(options, list):
            # List of methods - assume max_n_body = nfrags?
            # Empty list - assume name and max_n_body = nfrags?
            # List consisting of an int - assume int case
            raise Exception('List is not a currently implemented form of n_body '
                            'input. Try a dict.')
        else:
            func = db['n_body_func']
            for key in options.keys():
                # wfn method?
                if key in psi4.driver.procedures['{}'.format(func.__name__)].keys():
                    processed_options['methods'].update({key:options[key]})
                # dft method?
                elif (key == 'b3lyp'):
                    processed_options['methods'].update({key:options[key]})
#                # distributed?
#                elif key == 'distributed':
#                    print("distributed")
#                    processed_options['distributed'] = options[key]
#                elif key == 'cutoff':
#                    print("cutoff")
#                    processed_options['cutoff'] = options[key]
#                elif key == 'num_threads':
#                    print("num_threads")
#                    processed_options['num_threads'] = options[key]
                elif key == 'bsse':
#                    print("bsse")
                    processed_options['bsse'] = options[key]
#                    print('options[key]:')
#                    print(options[key])
#                elif key == 'pcm':
#                    print("pcm")
#                    processed_options['pcm'] = options[key]
                else:
                    raise Exception('Unrecognized n_body option {}.'.format(key))
                # remove the entry from options
#                del options[key]
#                print("key deleted")
    except Exception as E:
        print("Caught Exception E: {},{}".format(type(E)," ".join(E.args)))
        raise Exception('Unrecognized value for n_body input parameter. It may be '
                        ' a bool, int, list, or dict.')
            
    return processed_options


def cluster_dir(frags):
    directory = ''
    frags_copy = copy.deepcopy(frags)
    while frags_copy:
        directory += str(frags_copy.pop())
        if frags_copy:
            directory += '_'
    return(directory)


def n_body_dir(n):
    num_2_word = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven',
                  'eight', 'nine', 'ten', 'eleven', 'twelve', 'thirteen',
                  'fourteen', 'fifteen', 'sixteen', 'seventeen', 'eighteen',
                  'nineteen', 'twenty', 'twenty-one', 'twenty-two',
                  'twenty-three', 'twenty-four', 'twenty-five']

    if not isinstance(n, int):
        n = int(n[0])

    return('{}_body'.format(num_2_word[n]))

def harvest_data(db,method,n):
    for result in db[method]['results']:
        getattr(sys.modules[__name__],'harvest_{}_data'.format(result))(db,method,n)

def harvest_g09(db,method,n):
    print('Harvesting g09...')
    for result in db[method]['results']:
        print ("Harvesting {}".format(result))
        getattr(sys.modules[__name__],'harvest_g09_{}'.format(result))(db,method,n)

def harvest_scf_energy_data(db,method,n):
    body = n_body_dir(n)
    for job in db[method][n]['job_status']:
        with open('{}/{}/{}/output.json'.format(method,body,job),'r') as outfile:
            jout = json.load(outfile)
            energy = jout["SCF TOTAL ENERGY"]
            db[method][n]['scf_energy']['raw_data'].update({job:energy})

def harvest_scf_dipole_data(db,method,n):
    body = n_body_dir(n)
    for job in db[method][n]['job_status']:
        with open('{}/{}/{}/output.json'.format(method,body,job),'r') as outfile:
            jout = json.load(outfile)
            x = jout["SCF DIPOLE X"]
            y = jout["SCF DIPOLE Y"]
            z = jout["SCF DIPOLE Z"]
            db[method][n]['scf_dipole']['raw_data'].update({job:[x,y,z]})

def harvest_quadrupole_data(db,method,n):
    body = n_body_dir(n)
    for job in db[method][n]['job_status']:
        with open('{}/{}/{}/output.json'.format(method,body,job),'r') as outfile:
            jout = json.load(outfile)
            xx = jout["CC QUADRUPOLE XX"]
            xy = jout["CC QUADRUPOLE XY"]
            xz = jout["CC QUADRUPOLE XZ"]
            yy = jout["CC QUADRUPOLE YY"]
            yz = jout["CC QUADRUPOLE YZ"]
            zz = jout["CC QUADRUPOLE ZZ"]
            db[method][n]['quadrupole']['raw_data'].update({job:[xx,xy,xz,yy,yz,zz]})

def harvest_rotation_data(db,method,n):
    body = n_body_dir(n)
    name = method.upper()
    for omega in db['omega']:
        db[method][n]['rotation']['raw_data'][omega] = {}
    for job in db[method][n]['job_status']:
        with open('{}/{}/{}/output.json'.format(method,body,job),'r') as outfile:
            jout = json.load(outfile)
            # NOTE: Only worrying about Modified Velocity Gauge data for now
            for omega in db['omega']:
                optrot = jout["{} SPECIFIC ROTATION (MVG) @ {}NM".format(name, omega)]
                db[method][n]['rotation']['raw_data'][omega].update({job: optrot})
                

def harvest_rotation_tensor_data(db, method, n):
    body = n_body_dir(n)
    # Creates one massive list of all the elements from each tensor
    for job in db[method][n]['job_status']:
        tensors = []
        with open('{}/{}/{}/output.dat'.format(method,body,job),'r') as outfile:
            while True:
                # Get ONLY MVG data
                tensor = grab_psi4_matrix(outfile, 'Optical Rotation '
                                                'Tensor (Modified', 3)
                if not tensor:
                    break
                else:
                    tensors.extend(tensor)
        db[method][n]['rotation_tensor']['raw_data'].update({job: tensors})

def harvest_polarizability_data(db, method, n):
    body = n_body_dir(n)
    name = method.upper()
    for omega in db['omega']:
        db[method][n]['polarizability']['raw_data'][omega] = {}
    for job in db[method][n]['job_status']:
        with open('{}/{}/{}/output.json'.format(method,body,job),'r') as outfile:
            jout = json.load(outfile)
            # Need to fix psivars dict (called w/ core.get_variables(), 
            # stored in output.json) to include ALL omegas, not just 
            # the last one computed
            for omega in db['omega']:
                pols = jout["{} DIPOLE POLARIZABILITY @ {}NM".format(name, omega)]
                db[method][n]['polarizability']['raw_data'][omega].update({job: pols})

def harvest_polarizability_tensor_data(db, method, n):
    body = n_body_dir(n)
    # Creates one massive list of all the elements from each tensor
    for job in db[method][n]['job_status']:
        tensors = []
        with open('{}/{}/{}/output.dat'.format(method,body,job),'r') as outfile:
            while True:
                tensor = grab_psi4_matrix(outfile, 'Dipole Polarizability', 3)
                if not tensor:
                    break
                else:
                    tensors.extend(tensor)
        #Resort the tensors into user specified order 
        db[method][n]['polarizability_tensor']['raw_data'].update({job: tensors})

def harvest_cc_energy_data(db, method, n):
    body = n_body_dir(n)
    name = method.upper()
    for job in db[method][n]['job_status']:
        with open('{}/{}/{}/output.json'.format(method,body,job),'r') as outfile:
            jout = json.load(outfile)
            energy = jout["{} TOTAL ENERGY".format(name)]
            db[method][n]['cc_energy']['raw_data'].update({job:energy})

def harvest_cc_energy(db,method,n):
#####NOTE: is this function even useful??? find where this is used!!#####
    body = n_body_dir(n)
    name = psi4.get_global_option('WFN')
    for job in db[method][n]['job_status']:
        with open('{}/{}/{}/output.dat'.format(method,body,job)) as outfile:
            for line in outfile:
                if '{} correlation energy'.format(name.upper()) in line:
                    try:
                        (i,j,k,l, cc_energy) = line.split()
                        db[method][n]['cc_energy']['raw_data'].update({job: [float(cc_energy)]})
                    except:
                        pass

def harvest_g09_scf_energy(db,method,n):
    body = n_body_dir(n)
    for job in db[method][n]['job_status']:
        with open('{}/{}/{}/Test.FChk'.format(method,body,job)) as outfile:
            for line in outfile:
                if 'SCF Energy' in line:
                    (i, i, i, energy) = line.split()
                    db[method][n]['scf_energy']['raw_data'].update({job:float(energy)})

def harvest_g09_scf_dipole(db,method,n):
    body = n_body_dir(n)
    # FChk stores dipole moment in AU, convert to Debye
    factor = psi4.constants.dipmom_au2debye
    for job in db[method][n]['job_status']:
        get_next = False
        with open('{}/{}/{}/Test.FChk'.format(method,body,job)) as outfile:
            for line in outfile:
                if get_next:
                    (x, y, z) = line.split()
                    db[method][n]['scf_dipole']['raw_data'].update({job:[float(x)*factor,float(y)*factor,float(z)*factor]})
                    get_next = False
                if 'Dipole Moment' in line:
                    get_next = True


def harvest_g09_quadrupole(db,method,n):
    # ONLY HARVESTING THE TRACELESS QUADRUPOLE MOMENT!
    body = n_body_dir(n)
    # Convert AU to Debye-Angstrom
    factor = psi4.constants.dipmom_au2debye * psi4.constants.bohr2angstroms
    for job in db[method][n]['job_status']:
        get_next = 0
        with open('{}/{}/{}/Test.FChk'.format(method,body,job),'r') as outfile:
            for line in outfile:
                if get_next == 1:
                    (xx, yy, zz, xy, xz) = line.split()
                    get_next=2
                    continue
                if get_next == 2:
                    yz = line
                    db[method][n]['quadrupole']['raw_data'].update({job:[float(xx)*factor,float(yy)*factor,float(zz)*factor,float(xy)*factor,float(xz)*factor,float(yz)*factor]})
                    get_next = 0
                    
                if 'Quadrupole Moment' in line:
                    get_next = 1
#                if get_next == 1:
#                    (xx, xx_val, yy, yy_val, zz, zz_val) = line.split()
#                    get_next = 2
#                    continue
#                if get_next == 2:
#                    (xy, xy_val, xz, xz_val, yz_val, yz_val) = line.split()
#                    db[method][n]['quadrupole']['raw_data'].update({job: 
#                                [float(xx_val), float(yy_val), float(zz_val),
#                                 float(xy_val), float(xz_val), float(yz_val)]})
#                    get_next = 0
#                if 'Traceless Quadrupole moment' in line:
#                    get_next = 1

def harvest_g09_rotation(db,method,n):
    """Harvests g09 rotations. I have to extend the database manually here, so that
        harvest_g09_rotation_tensor() doesn't get called multiple times."""
    db[method]['results'].append('rotation_tensor')
    db[method][n]['rotation_tensor'] = collections.OrderedDict()
    db[method][n]['rotation_tensor']['correction'] = 0
    db[method][n]['rotation_tensor']['rotation_tensor'] = 0
    db[method][n]['rotation_tensor']['vmfc_correction'] = 0
    db[method][n]['rotation_tensor']['vmfc_approximation'] = 0
    db[method][n]['rotation_tensor']['mbcp_correction'] = 0
    db[method][n]['rotation_tensor']['mbcp_approximation'] = 0
    db[method][n]['rotation_tensor']['raw_data'] = collections.OrderedDict()
    db[method][n]['rotation_tensor']['cooked_data'] = collections.OrderedDict()
    harvest_g09_rotation_tensor(db,method,n)
    body = n_body_dir(n)
    c = psi4.constants.c
    h = psi4.constants.h
    h2j = psi4.constants.hartree2J
    Na = psi4.constants.na
    me = psi4.constants.me
    hbar = h / 2.0 / math.pi
    prefactor = -72E6 * hbar**2 * Na / c**2 / me**2

    for job in db[method][n]['job_status']:
        M = db[method][n]['MW'][job]
        optrot = []
        for omega in db[method][n]['rotation_tensor']['raw_data'][job]:
            w_h = c * h * 1E9 / h2j / omega
            tr = np.trace(db[method][n]['rotation_tensor']['raw_data'][job][omega])
            # Multiply by additional omega to account for different mu operator
            rot = prefactor * w_h * tr / M / 3.0 * w_h
            optrot.append(rot)
        optrot = reorder_g09_rotations(optrot, db)
        db[method][n]['rotation']['raw_data'].update({job: optrot})
        
#        with open('{}/{}/{}/input.log'.format(method,body,job)) as outfile:
#            for line in outfile:
#                # Assume specific rotations are presented in descending order in
#                # terms of wavelength in nm, no matter the user specified
#                # order.
#                if '[Alpha]' in line:
#                    try:
#                        (i,i,i,i,i,i,i,i,i,i,val,i) = line.split()
#                        optrot.append(float(val))
#                    # if values are too large there is no space between = and result
#                    except ValueError:
#                        try:
#                            (i,i,i,i,i,i,i,i,i,val,i) = line.split()
#                            optrot.append(float(val[1:]))
#                        except:
#                            optrot.append(0.00)
#                            print('There has been an overflow in the optical rotation data and they are now meaningless after {}-body.'.format(n))
        # Resort the rotations from descending wavelength (nm) order to user
        # specified order
#        optrot = reorder_g09_rotations(optrot, db)

        # Add list of rotations onto database entry
#        print("Got the rotations! We found: {} for job {}".format(optrot, job))
#        db[method][n]['rotation']['raw_data'].update({job: optrot})

def reorder_g09_rotations(optrot, db, omega=None):
    '''# g09 specific rotations are output in order of wavelength
    # ascending/descending depending on the units
    # Descending:  NM
    # Ascending: AU, EV, HZ 
    # No matter the user specified order. Here we're resorting them to
    # be in the same order as specified by the user.'''

    # Get user specified omega list
#    if not omega:
#        omega = psi4.get_global_option('OMEGA')

#    if len(omega) > 1:
#        # Remove the units entry from omega
#        units = omega.pop()
#    else:
#        units = 'au' #default units
    # just working with units in nm for now
    # last entry is always the units
    if 'omega' in db:
        omega_list = copy.deepcopy(db['omega'])
#        units = omega_list.pop()

    
    # Omega order contains a list of numbers for keeping track of the
    # original user specified order of omega
    omega_order = []
    for n in range(1,len(omega_list)+1):
        omega_order.append(n)

    # Zip omega with our order tracking list
    # sort the zipped list in terms of omega (descending, assuming NM)
    temp = sorted(zip(omega_list, omega_order),reverse=True)

    # Unzip the sorted lists @omega_descending contains the wavelengths in
    # descending order @omega_order contains the reordered list of numbers
    # allowing us to get back to user entered order
    omega_descend, omega_order = zip(*temp)

    # Check if the data is tensors (i.e. longer than the number of omegas)
    if len(optrot) > len(omega_list):
        temp_data = []
        for i in range(len(omega_list)):
            tensor = []
            # Assuming 3x3 tensors
            for j in range(9):
                tensor.append(optrot.pop(0))
            temp_data.append(tensor)
        optrot = temp_data
    # Zip together our list of the original order and the collected rotations
    # Sort the zipped list according to omega_order
    temp = sorted(zip(omega_order, optrot))
    # Unzip to get the rotations in the user specified order
    junk, optrot = zip(*temp)
    if isinstance(optrot[0],list):
        temp_data = []
        for n in optrot:
            temp_data.extend(n)
        optrot = temp_data
    else:
        # convert from tuple to list
        optrot = [n for n in optrot]
    # Return the list
    return optrot


def harvest_g09_rotation_tensor(db, method, n):
    body = n_body_dir(n)
    # G09 ALWAYS prints rotation for each wavelength in DESCENDING order
    omega_dec = copy.deepcopy(db['omega'])
    omega_dec.sort(reverse=True)
    n_omega = len(omega_dec)
    # Gaussian prints every single value from each tensor in a row, 5 values per row
    n_rows = math.floor(9*n_omega/5)
    for job in db[method][n]['job_status']:
        # Need a dictionary to hold the different tensors
        tensors = {}       
        # Need to cram all of the values into a 1D list
        ten_vals = []
        get_line = 0
        get_next = 0
        with open('{}/{}/{}/Test.FChk'.format(method,body,job)) as outfile:
            for line in outfile:
                if (get_next == 1) & (get_line <= n_rows):
                    ten_vals += line.split()
                    get_line += 1
                if 'FD Optical Rotation Tensor' in line:
                    get_next = 1
        
        for i in range(0,len(ten_vals)):
            ten_vals[i] = float(ten_vals[i])
        vals = np.asarray(ten_vals)
        vals = np.split(vals, n_omega)
        i = 0
        for omega in omega_dec:
            tensors[omega] = vals[i].reshape(3,3).T
            i+=1
#        # Resort the tensors to user specified order
#        tensors = reorder_g09_rotations(tensors)
        db[method][n]['rotation_tensor']['raw_data'].update({job: tensors})
            

def harvest_g09_polarizability(db, method, n):
    body = n_body_dir(n)
    omega = psi4.get_global_option('OMEGA')
    if len(omega) > 1:
        omega.pop()
    n_omega = len(omega)
    for job in db[method][n]['job_status']:
        pols = []
        with open('{}/{}/{}/input.log'.format(method,body,job),'r') as outfile:
            for line in outfile:
                # polarizabilities are presented in order of increasing omega
                # (in au). Means I should probably pass them through the sorting
                # function I use for rotations
                if 'Isotropic polarizability' in line:
                    (i,i,i,i,i, value, i) = line.split()
                    pols.append(float(value))
        # Remove zero frequency result
        if len(pols) > n_omega:
            pols.pop(0)    

        # Resort the rotations from descending wavelength (nm) order to user
        # specified order
        pols = reorder_g09_rotations(pols)

        # Place results in database
        db[method][n]['polarizability']['raw_data'].update({job: pols})


                
def harvest_g09_polarizability_tensor(db, method, n):
    body = n_body_dir(n)
    omega = psi4.get_global_option('OMEGA')
    if len(omega) > 1:
        omega.pop()
    n_omega = len(omega)
    for job in db[method][n]['job_status']:
        tensors = []
        with open('{}/{}/{}/input.log'.format(method,body,job),'r') as outfile:
            for k in range(1, n_omega+1):
                # n+1 skips the zero tensor which is _frequency 1_
                tensors.extend(grab_g09_matrix(outfile, 'Alpha(-w,w) frequency  {}'.format(k+1), 3))
        tensors = reorder_g09_rotations(tensors)
        db[method][n]['polarizability_tensor']['raw_data'].update({job: tensors})

def cook_data(db, method, n):
    # Automatic cooking requires all result data be stored as lists in the
    # database. Even if the result is a single number (i.e. energies)
    cooked_data = collections.OrderedDict()
    raw_data    = collections.OrderedDict()
    for result in db[method]['results']:
        # Read in all cooked_data for m < n
        for m in range(1, n):
            cooked_data[m] = db[method][m][result]['cooked_data']

        # Read in raw data for m
        # Eventually should just read this in to cooked_data[n]
        #cooked_data[n] = db[method][n][result]['raw_data']
        raw_data[n] = db[method][n][result]['raw_data']

        # Cook n data
        cooked_data[n] = copy.deepcopy(raw_data[n])
        for m in range(1,n):
            print(result)
            print(cooked_data[m])
            for n_key in cooked_data[n]:
                for m_key in cooked_data[m]:
                    m_set = set(m_key.split('_'))
                    n_set  = set(n_key.split('_'))
                    if len(m_set.intersection(n_set)) == len(m_set):
                        cooked_data[n][n_key] = [x-y for x,y in 
                        zip(cooked_data[n][n_key],cooked_data[m][m_key])]

        # Set correction list length
        #try:
        #    correction_length = len(cooked_data[n].itervalues().next())
        #except TypeError:
            
        # There is a problem with this when the data doesn't actually exist
        # (i.e. when calling property('scf', properties=[polarizability] because
        # there is no polarizability so the list is NoneType. The above try fix
        # might work, but it would be better to figure out how to actually know
        # what data quantities are available in a smart way.
#        correction = [0.0 for i in
#                      range(len(cooked_data[n].itervalues().next()))]
        correction = [0.0 for i in
                      range(len(iter(cooked_data[n].values()).next()))]
        # Add up all contributions to current correction
        for key,val in cooked_data[n].items():
            correction = [x+y for x,y in zip(correction,val)]

        db[method][n][result]['cooked_data'] = cooked_data[n]
        db[method][n][result]['correction']  = copy.deepcopy(correction)
        # Final result value is the correction for n plus the result value from
        # n-1 
        db[method][n][result][result] = copy.deepcopy(correction)
        if n > 1: 
            prev = db[method][n-1][result][result]
            curr = db[method][n][result][result]
            db[method][n][result][result] = [x+y for x,y in zip(curr,prev)]

        # Wipe out data dicts for next result
        cooked_data.clear()
        raw_data.clear()

def mbcp_cook(db, method, n):
    cooked_data = collections.OrderedDict()
    raw_data = collections.OrderedDict()

    for result in db[method]['results']:
        # Create list of required cooked data
        garden = [1]
        print("MBCP cook field = {}".format(n))
        #for n in range(2, field):
        for m in range(2, n):
            garden.append('{}-{}r'.format(m,1))

        field = '{}-{}r'.format(n,1)
        print(field)
        print(result)

        for bed in garden:
            cooked_data[bed] = db[method][bed][result]['cooked_data']

        cooked_data[field] = copy.deepcopy(db[method][field][result]['raw_data'])

        for bed in garden:
            for field_key in cooked_data[field]:
                for bed_key in cooked_data[bed]:
                    field_set = set(field_key.split('_'))
                    bed_set = set(bed_key.split('_'))
                    if len(bed_set.intersection(field_set)) == len(bed_set):
                        cooked_data[field][field_key] = [y-x for x,y in
                                          zip(cooked_data[field][field_key],cooked_data[bed][bed_key])]
        correction = [0.0 for i in
                      range(len(cooked_data[field].itervalues().next()))]

        for key,val in cooked_data[field].items():
            correction = [x+y for x,y in zip(correction,val)]

        # Individual counterpoise corrections
        db[method][field][result]['cooked_data'] = cooked_data[field]
        # Total counterpoise correction
        db[method][n][result]['mbcp_correction']  = copy.deepcopy(correction)

        # This contains the counterpoise corrected approximation i.e. naive
        # n-body result minus the MBCP(n) correction
        db[method][n][result]['mbcp_approximation'] = copy.deepcopy(correction)
        db[method][n][result]['mbcp_approximation'] = [x+y for x,y in
                zip(db[method][n][result][result],db[method][n][result]['mbcp_approximation'])]

        # Wipe out data dicts for next result
        cooked_data.clear()
        raw_data.clear()

def vmfc_cook(db, method, n):
    # Automatic cooking requires all result data be stored as lists in the
    # database. Even if the result is a single number (i.e. energies)
    cooked_data = collections.OrderedDict()
    raw_data    = collections.OrderedDict()

    for result in db[method]['results']:

        # Read in all levels m < n in the n-mer basis
        for m in range(1,n):
            cooked_data['{}-{}r'.format(n,m)] = db[method]['{}-{}r'.format(n,m)][result]['raw_data']
        # Read in n-mers in n-mer basis
        cooked_data[n] = db[method][n][result]['raw_data']

        # Preparation of all V's less than Vn
        for v in range(2,n):
            for m in range(1,v):
                for n_key in cooked_data['{}-{}r'.format(n,v)]:
                    for m_key in cooked_data['{}-{}r'.format(n,m)]:
                        m_set = set(m_key.split('_'))
                        n_set = set(n_key.split('_'))
                        mg_set = set(re.sub('g','',m_key).split('_'))
                        ng_set = set(re.sub('g','',n_key).split('_'))
                        if len(m_set.intersection(n_set)) == len(m_set)-(v-m) and len(mg_set.intersection(ng_set)) == n:
                            cooked_data['{}-{}r'.format(n,v)][n_key] = [x-y for x,y in 
                            zip(cooked_data['{}-{}r'.format(n,v)][n_key],cooked_data['{}-{}r'.format(n,m)][m_key])]


        # Cook n data
        for m in range(1,n):
            for n_key in cooked_data[n]:
                for m_key in cooked_data['{}-{}r'.format(n,m)]:
                    m_set  = set(re.sub('g','',m_key).split('_'))
                    n_set = set(n_key.split('_'))
                    if len(m_set.intersection(n_set)) == n:
                        cooked_data[n][n_key] = [x-y for x,y in 
                        zip(cooked_data[n][n_key],cooked_data['{}-{}r'.format(n,m)][m_key])]

        # Set correction list length
        #try:
        #    correction_length = len(cooked_data[n].itervalues().next())
        #except TypeError:
            
        # There is a problem with this when the data doesn't actually exist
        # (i.e. when calling property('scf', properties=[polarizability] because
        # there is no polarizability so the list is NoneType. The above try fix
        # might work, but it would be better to figure out how to actually know
        # what data quantities are available in a smart way.
        correction = [0.0 for i in
                      range(len(cooked_data[n].itervalues().next()))]
        # Add up all contributions to current correction
        for key,val in cooked_data[n].items():
            correction = [x+y for x,y in zip(correction,val)]

        db[method][n][result]['vmfc_correction']  = copy.deepcopy(correction)
        # Final result value is the correction for n plus the result value from
        # n-1 
        db[method][n][result]['vmfc_approximation'] = copy.deepcopy(correction)
        if n > 1: 
            prev = db[method][n-1][result]['vmfc_approximation']
            curr = db[method][n][result]['vmfc_approximation']
            db[method][n][result]['vmfc_approximation'] = [x+y for x,y in zip(curr,prev)]

        # Wipe out data dicts for next result
        cooked_data.clear()
        raw_data.clear()


def grab_psi4_matrix(outfile, matrix_name, row_tot):
    collect_matrix = False
    n_rows = 0
    n_tries = 0
    matrix_data = []
    for line in outfile:
        if matrix_name in line:
            collect_matrix = True
        if collect_matrix and (n_rows < row_tot):
            try:
                n_tries += 1
                if n_tries > (row_tot + 13):
                    raise Exception('{} matrix was unreadable. Scanned {} '
                                    'lines'.format(matrix_name, n_tries))
                else:
                    (index, x, y, z) = line.split()
                    matrix_data.append(float(x))
                    matrix_data.append(float(y))
                    matrix_data.append(float(z))
                    n_rows += 1
            except:
                pass

        if (n_rows == row_tot) and (len(matrix_data) != 3*row_tot):
            raise Exception('Collecting matrix data failed!')

        if len(matrix_data) == 3*row_tot:
            return matrix_data
    else:
        return None

def grab_g09_matrix(outfile, matrix_name, row_tot):
    collect_matrix = False
    n_rows = 0
    n_tries = 0
    matrix_data = []
    for line in outfile:
        if matrix_name in line:
            collect_matrix = True
        if collect_matrix and (n_rows < row_tot):
            try:
                n_tries += 1
                if n_tries > (row_tot + 5):
                    raise Exception('{} matrix was unreadable. Scanned {} '
                                    'lines'.format(matrix_name, n_tries))
                else:
                    (index, x, y, z) = line.split()
                    # Change g09 exponent style to python style
                    matrix_data.append(float(x.replace('D','E')))
                    matrix_data.append(float(y.replace('D','E')))
                    matrix_data.append(float(z.replace('D','E')))
                    n_rows += 1
            except:
                pass

        if (n_rows == row_tot) and (len(matrix_data) != 3*row_tot):
            raise Exception('Collecting matrix data failed!')

        if len(matrix_data) == 3*row_tot:
            return matrix_data
    else:
        return None

def interacting_pairs(cutoff):
    nfrag = 0
    frag_list = [[]]
    i_pairs = set()
    with open('input.dat') as infile:
        for line in infile:
            try:
                (atom, x, y, z) = line.split()
                frag_list[nfrag].append([float(x), float(y), float(z)])
            except:
                if '--' in line:
                    # Creates the entry in frag_list for next fragments coordinates
                    frag_list.append([])
                    # increment frag counter
                    nfrag += 1

    for n,frag in enumerate(frag_list):
        for n2,frag2 in enumerate(frag_list[n+1:]):
            for atom in frag:
                for atom2 in frag2:
                    dist = distance(atom, atom2)
                    if dist < cutoff:
                        #i_pairs.add((n+1, n2+1+(n+1)))
                        i_pairs.add((n2+1+(n+1), n+1))
    #return sorted(list(i_pairs),reverse=True)
    return sorted([list(pair) for pair in i_pairs])


def distance(p1, p2):
    x_diff = p1[0] - p2[0]
    y_diff = p1[1] - p2[1]
    z_diff = p1[2] - p2[2]
    return math.sqrt(x_diff**2 + y_diff**2 + z_diff**2)

# Include atoms within basis_cutoff of cluster as ghosts
def expand_basis(clusters, indexes, basis_cutoff):
    # Set up fragment coordinates
    nfrag = 0
    frag_list = [[]]
    with open('input.dat') as infile:
        for line in infile:
            try:
                (atom, x, y, z) = line.split()
                frag_list[nfrag].append([float(x), float(y), float(z)])
            except:
                if '--' in line:
                    # Creates the entry in frag_list for next fragments coordinates
                    frag_list.append([])
                    # increment frag counter
                    nfrag += 1

    for cluster in clusters:
        print(cluster)
        #for atom in cluster

        

def cluster_distance(n,index):

    # Set up fragment coordinates
    nfrag = 0
    frag_list = [[]]
    with open('input.dat') as infile:
        for line in infile:
            try:
                (atom, x, y, z) = line.split()
                frag_list[nfrag].append([float(x), float(y), float(z)])
            except:
                if '--' in line:
                    # Creates the entry in frag_list for next fragments coordinates
                    frag_list.append([])
                    # increment frag counter
                    nfrag += 1

    # Compute interfragment distances
    if len(index) == 2:
        m = index.pop()
        k = index.pop()
        min_dist = 1000.0
        for atom in frag_list[m-1]:
            for atom2 in frag_list[k-1]:
                dist = distance(atom,atom2)
                if dist < min_dist:
                    min_dist = dist

        return min_dist

def plant(cluster, db, kwargs, method, directory):
    # Create job directory
    try:
        os.makedirs(directory)
    except:
        if os.path.isdir(directory):
            print(directory)
            pass
        else:
            raise Exception('Attempt to create {} directory '
                            'failed.'.format(directory))
    # Write psi4 input files
    # What type of function are we running?
    func = db['n_body_func']
    if (method == 'b3lyp'):
        # Using g09 for dft properties
        # Convert from psi4 keywords to g09 equivalent
        # Note that Polar=Dipole calculates quadrupoles as well
        psi4_to_g09 = { 'rotation': 'Polar=OptRot',
                        'polarizability': 'Polar=Dipole',
                        'quadrupole': 'Polar=Dipole' }
        # Did the user specify omegas?
#        omega_list = psi4.get_global_option('OMEGA')
        if 'omega' in db:
            omega_list = db['omega']
            units = 'nm'
        # NOTE: doing everything in nm for now... should fix this
#        if len(omega_list) > 1:
#            units = omega_list.pop()
#        else:
#            units = 'AU'

        infile = open('{}/input.com'.format(directory),'w')
        basis = psi4.core.get_global_option('BASIS')
#        basis = psi4.core.print_option('BASIS')
        # Check for n-aug-cc-pVXZ (G09 doesn't have them)
        # GEN basis is hard wired for n-aug-cc-pVXZ
        real_basis = ''
        gen_basis_list = ['D-AUG-CC-PVDZ','T-AUG-CC-PVDZ','Q-AUG-CC-PVDZ',
                          'D-AUG-CC-PVTZ','T-AUG-CC-PVTZ','Q-AUG-CC-PVTZ',
                          'D-AUG-CC-PVQZ','T-AUG-CC-PVQZ','Q-AUG-CC-PVQZ']
        if basis in gen_basis_list:
            real_basis = copy.deepcopy(basis)
            basis = 'GEN'
        # link keeps up with running multiple jobs in single
        # input.com file
        link = False
        # Write the job line command
        if 'properties' in kwargs:
            for prop in kwargs['properties']:
                if link:
                    infile.write('\n--Link1--\n\n')
                if db['num_threads']:
                    infile.write('%NProcShared={}\n'.format(db['num_threads']))
                infile.write('%Mem={}MB\n'.format(int(psi4.get_memory()/1000000)))
                # G09 only likes integer values on MEM
                infile.write('#p {}/{} NoSymmetry'.format(method,basis))
                if omega_list:
                    infile.write(' CPHF(RdFreq)')
                    # Uncomment following line for tight convergence
                    #infile.write(' CPHF(Conver=12,RdFreq,Grid=UltraFine)')
                if db['pcm']:
                    infile.write(' SCRF=(PCM,Solvent={})'.format(db['pcm']))
                # Uncomment following line for tight convergence
                #infile.write(' SCF(Conver=12,MaxCycle=512) Integral(UltraFine)')
                infile.write(' {} FormCheck\n\n'.format(psi4_to_g09[prop]))

                # Write autogen comment and job name info
                infile.write('This is a g09 input file auto-generated '
                         'by psi4 for an n_body {} job.\n'.format(func.__name__))
                infile.write('{}\n\n'.format(cluster.name()))

                # Write geometry (units = ang)
#                infile.write('{}\n'.format(cluster.save_string_xyz_g09()))
#                infile.write('{}\n'.format(cluster.save_string_xyz(False)))
                infile.write('{}\n'.format(save_geom_string_safe(cluster,'g09')))
                # Write user provided wavelengths
                for wavelength in omega_list:
                    if units.upper() == 'NM':
                        infile.write('{}nm '.format(wavelength))
                    elif units.upper() == 'AU':
                        infile.write('{}au '.format(wavelength))
                    else:
                        # Would be good to auto-convert from hz and ev
                        # instead, but until then
                        raise ValidationError('Wavelengths must be in '      
                                              'NM or AU for use with g09')
                infile.write('\n')
                # Check for GEN basis
                # Works only with n-aug-cc-pVXZ basis
                if basis == 'GEN':
                    infile.write('\n@{}/basis_sets/gaussian/{}.gbs'.format(os.getenv('HOME'),real_basis))
                link = True
        # Otherwise assume energy calculation
        else:
            if db['num_threads']:
                infile.write('%NProcShared={}\n'.format(db['num_threads']))
            infile.write('%Mem={}MB\n'.format(int(psi4.get_memory()/1000000)))
            infile.write('#p {}/{} NoSymmetry FormCheck'.format(method,basis))
            if db['pcm']:
                infile.write(' SCRF=(PCM,Solvent={})'.format(db['pcm']))
            # Uncomment following line for tight convergence
            #infile.write(' SCF(Conver=12,MaxCycle=512) Integral(UltraFine)')

            # Write autogen comment and job name info
            infile.write('\n\nThis is a g09 input file auto-generated '
                     'by psi4 for an n_body {} job.\n'.format(func.__name__))
            infile.write('{}\n\n'.format(cluster.name()))

            # Write geometry (units = ang)
#            infile.write('{}\n'.format(cluster.save_string_xyz_g09()))
            #infile.write('{}\n'.format(cluster.save_string_xyz(False)))
#            print(cluster)
#            print("Here is the output: {}".format(save_geom_string_safe(cluster,'g09')))
            infile.write('{}\n'.format(save_geom_string_safe(cluster,'g09')))
            # Write user provided wavelengths
            infile.write('\n')
            # Check for GEN basis
            # Hardwired for daDZ for now
            if basis == 'GEN':
                infile.write('\n@{}/basis_sets/gaussian/daDZ.gbs'.format(os.getenv('HOME')))

        infile.write('\n')
        infile.close()
#    # Otherwise assume method is a wfn method. All methods were
#    # checked earlier for validity.
    else:
        infile = open('{}/input.dat'.format(directory),'w')
        infile.write('# This is a psi4 input file auto-generated for an'
                 ' n_body {} job.\n'.format(func.__name__))
        mol_open = 'molecule ' + cluster.name() + ' {\n'
        mol_close = '}'
        infile.write("{}{}{}".format(mol_open,cluster.create_psi4_string_from_molecule(),mol_close))
        infile.write(psi4.p4util.format_options_for_input())
        infile.write("\n{}('{}', ".format(func.__name__, method))
        for key,val in kwargs.items():
        #infile.write('{}={}, '.format(key,val))
            infile.write('{}={} '.format(key,val))
        infile.write(')\n')
        infile.write("""\n\nimport json\nwith open('output.json', 'w') as dumpf:\n""")
        infile.write("""    json.dump(core.get_variables(), dumpf, indent=4)""")
        infile.close()

def ghost_dir(indexes, ghost):
    directory = ''
    indexes_copy = copy.deepcopy(indexes)
    while indexes_copy:
        index = indexes_copy.pop()
        directory += str(index)
        if index in list(ghost):
            directory += 'g'
        if indexes_copy:
            directory += '_'
    return(directory)

def save_geom_string_safe(mol, mode):
    """Returns the cartesian geometry of a molecule as a string.
      Parameters
      ----------
      mol : psi4.core.Molecule object
      mode : string, print mode either 'psi4' or 'g09'

      Returns
      -------
      geom : string, formatted cartesian geometry suitable for input file(s) 

      Notes
      -----
      Please use Angstroms for now. I'm not sure how other units will react.
    """

    # mol.fn(i) gives in bohr, so must convert back to angstroms if input file asks
    factor = psi4.constants.bohr2angstroms if str(mol.units) == 'GeometryUnits.Angstrom' else 1.0 

    geom = '{},{}\n'.format(mol.molecular_charge(),mol.multiplicity())

    if mode.lower() == 'psi4':
        for i in range(mol.natom()):
            geom += ' {} {} {} {}\n'.format((mol.symbol(i) if mol.Z(i) else 'Gh({})'.format(mol.symbol(i))), mol.fx(i) * factor, mol.fy(i) * factor, mol.fz(i) * factor) 
            

    if mode.lower() == 'g09':
        for i in range(mol.natom()):
            geom += ' {} {} {} {}\n'.format((mol.symbol(i) if mol.Z(i) else '{}-Bq'.format(mol.symbol(i))), mol.fx(i) * factor, mol.fy(i) * factor, mol.fz(i) * factor) 

    return geom
