#
# @BEGIN LICENSE
#
# n_body by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util

##### IMPORTS REQUIRED FOR TAYLOR'S run_n_body() FUNCTION (SEE OLD proc.py)#####
import shelve
from . import n_body
import itertools
#os is already imported in n_body
import os

def run_n_body(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    n_body can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('n_body')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

#    print("START run_n_body()")

    # Your plugin's psi4 run sequence goes here
    psi4.core.set_local_option('MYPLUGIN', 'PRINT', 1)

    ##### BEGIN TAYLOR'S run_n_body() FUNCTION #####

    # Only energy and properties are implemented for now, not sure others even
    # make sense anyways.

    db = shelve.open('database',writeback=True)
#    print(db)
    # Initialize database
    if not 'initialized' in db:
        n_body.initialize_database(db, kwargs)
        db['initialized'] = True
    elif not db['initialized']:
        n_body.initialize_database(db, kwargs)
        db['initialized'] = True

    # Remove n_body_func, it's being tracked in the database and we don't want to propagate it
    if 'n_body_func' in kwargs:
        kwargs.pop('n_body_func')

    # Process user requested options
    n_body_options = n_body.process_options(name, db, kwargs.pop('n_body'))

    # Compare database options to n_body_options
    if 'distributed' in n_body_options:
        db['distributed'] = n_body_options['distributed']
        if not db['distributed']:
            raise Exception("Currently n_body jobs must be run in distributed mode"
                            " use the n_body wrapper instead.\n")
    if 'cutoff' in n_body_options:
        db['cutoff'] = n_body_options['cutoff']
    if 'num_threads' in n_body_options:
        db['num_threads'] = n_body_options['num_threads']
    if 'bsse' in n_body_options:
        db['bsse'] = n_body_options['bsse']
    if 'pcm' in n_body_options:
        db['pcm'] = n_body_options['pcm']
    # methods consistency check
    if 'methods' in n_body_options:
        for key,val in n_body_options['methods'].items():
            if key in db['methods'].keys():
                # requested n_body_max smaller than previously requested
                if val < db['methods'][key]:
                    raise Exception('n_body_max for {} has '
                                    'decreased.'.format(key))
                # requested n_body_max has increased
                elif val > db['methods'][key]:
                    # save requested n_body_max
                    db['methods'][key] = val
                    # update database entries
                    n_body.extend_database(db,kwargs)
                    # set flag to regenerate inputs
                    db['inputs_generated'] = False
            else: # method not already in database
                # add method and n_body_max to database
                db['methods'][key] = val
                # update database entries
                n_body.extend_database(db,kwargs)
                # set flag to regenerate inputs
                db['inputs_generated'] = False

    # Get complete system
    molecule = psi4.core.get_active_molecule()

    # Determine interacting fragments
    if db['cutoff']:
        i_pairs = n_body.interacting_pairs(db['cutoff'])
        i_clusters = []

    if not db['inputs_generated']:
        # Create method directories
        for method,n_body_max in db['methods'].items():
            try:
                os.makedirs(method)
            except:
                if os.path.isdir(method):
                    print(method)
                    pass
                else:
                    raise Exception('Attempt to create {} directory '
                                    'failed.'.format(method))
            # Create n_body_level directories
            for n in range(1, n_body_max+1):
                directory = '{}/{}'.format(method,n_body.n_body_dir(n))
                try:
                    os.makedirs(directory)
                except:
                    if os.path.isdir(directory):
                        print(directory)
                        pass
                    else:
                        raise Exception('Attempt to create {} directory '
                                        'failed.'.format(directory))
                # Create job input files
                
                # Grab all ghost functions if SSFC 
                # Else, get all possible clusters without ghost atoms
                if db['bsse'] == 'ssfc':
                    clusters = psi4.extract_clusters(molecule, True, n)
                else:
                    clusters = psi4.extract_clusters(molecule, False, n)

                # Flip distance-dropoff off for now
                #if db['bsse'] == 'radius'
                #clusters = extract_clusters(molecule, False, n)
                #clusters = n_body.expand_basis(clusters, db['basis_cutoff'])

                # Get list of fragments involved in clusters
                indexes = psi4.extract_cluster_indexing(molecule, n)
                print('Length of index array = {}'.format(len(indexes)))

                # Testing for extension to inclusion of nearby atoms (ghosted)
                #cluster_string = clusters[0].save_string_xyz_g09()
                #print(cluster_string)
                #cluster_string += molecule.extract_subsets(2).save_string_xyz_g09()
                #print(molecule.extract_subsets(2).save_string_xyz_g09())
                #print(clusters[0].save_string_xyz_g09())
                #print(cluster_string)
                #cluster_string = ''

                # Get interacting clusters
                if db['cutoff']:
                    if n == 1:
                        # Generate all inputs
                        pass
                    elif n == 2:
                        # generate inputs for i_pairs only
                        for k in range(len(clusters)-1, -1, -1):
                            clus = indexes[k]
                            if clus in i_pairs:
                                i_clusters.append(clus)
                            else:
                                # Remove non-interacting clusters from lists
                                clusters.pop(k)
                                indexes.pop(k)

                    elif n > 2:
                        prev_clusters = copy.deepcopy(i_clusters)
                        i_clusters = []

                        for k in range(len(clusters)-1, -1, -1):
                            # Check if interacting cluster
                            clus = indexes[k]
                            for p_clus in prev_clusters:
                                # previous cluster is contained within current cluster
                                if set(p_clus).issubset(set(clus)):
                                    # Get remaining frag number
                                    if len(set(clus).difference(set(p_clus))) == 1:
                                        l = set(clus).difference(set(p_clus)).pop()
                                    else:
                                        raise Exception('Length of difference '
                                                        'set is greater than 1')
                                    for m in p_clus:
                                        # Is this an interacting pair?
                                        if sorted([l, m], reverse=True) in i_pairs:
                                            i_clusters.append(clus)
                                            # Don't check anymore pairs
                                            break
                                    if clus in i_clusters:
                                        break
                            else:
                                # Remove non-interacting clusters from lists
                                clusters.pop(k)
                                indexes.pop(k)

                # Create input files
                db[method][n]['total_num_jobs'] = len(clusters)
                print('n = {}'.format(n))
                print('Number of jobs created = {}'.format(len(clusters)))
                for k in range(len(clusters)):
                    psi4.activate(clusters[k])
                    cluster = psi4.core.get_active_molecule()
#                    db[method][n]['MW'][job] = value
                    ### Set-up things for input generation
                    # Get list of combinations to ghost
                    # MBCP 
                    if db['bsse'] == 'mbcp' and n > 1:
                        mbcp = list(itertools.combinations(indexes[k],n-1))
                        for ghost in mbcp:
                            directory = '{}/{}/{}'.format(method,
                                                          n_body.n_body_dir(n),n_body.ghost_dir(indexes[k],ghost))
                            cluster.set_ghost_fragments(list(ghost))
                            cluster.update_geometry()
                            cluster.set_name('{}_{}'.format(molecule.name(),n_body.ghost_dir(indexes[k],ghost)))
                            n_body.plant(cluster, db, kwargs, method, directory)
                            # Update job_status dict
                            n_basis = n
                            n_real = n - len(ghost)
                            ghost_jobs = '{}-{}r'.format(n_basis,n_real)
#                            db[method][ghost_jobs]['total_num_jobs'] = len(mbcp) * len(clusters)
#                            db[method][ghost_jobs]['job_status'].update({n_body.ghost_dir(indexes[k],ghost):
#                                                                        'not_started'})
                            db[method][n]['total_num_jobs'] = len(mbcp) * len(clusters)
                            db[method][n]['job_status'].update({n_body.ghost_dir(indexes[k],ghost):
                                                                        'not_started'})
                            # Calculate molecular weight
                            totalmass = sum(mol2.mass(i) for i in range(mol2.natom()) if mol2.Z(i))
                            totmass = sum(cluster.mass(i) for i in range(cluster.natom()) if cluster.Z(i))
                            db[method][n]['MW'].update({n_body.ghost_dir(indexes[k],ghost):totmass})
                            # Reactivate fragments for next round
                            cluster.set_active_fragments(list(ghost))
                            cluster.update_geometry()
                    # VMFC
                    elif db['bsse'] == 'vmfc':
                        vmfc = []
                        for n_ghost in range(n-1, 0, -1):
                            vmfc.extend(list(itertools.combinations(indexes[k],n_ghost)))
                        for ghost in vmfc:
                            directory = '{}/{}/{}'.format(method,
                                                          n_body.n_body_dir(n),n_body.ghost_dir(indexes[k],ghost))
                            cluster.set_ghost_fragments(list(ghost))
                            cluster.update_geometry()
                            cluster.set_name('{}_{}'.format(molecule.name(),n_body.ghost_dir(indexes[k],ghost)))
                            n_body.plant(cluster, db, kwargs, method, directory)
                            # Update job_status dict
                            n_basis = n
                            n_real = n - len(ghost)
                            ghost_jobs = '{}-{}r'.format(n_basis,n_real)
                            #db[method][ghost_jobs]['total_num_jobs'] = len(vmfc) * len(clusters)
#                            db[method][ghost_jobs]['total_num_jobs'] = len(list(itertools.combinations(range(0,n),n_real))) * len(clusters)
#                            db[method][ghost_jobs]['job_status'].update({n_body.ghost_dir(indexes[k],ghost):
#                                                                        'not_started'})
                            db[method][n]['total_num_jobs'] = len(list(itertools.combinations(range(0,n),n_real))) * len(clusters)
                            db[method][n]['job_status'].update({n_body.ghost_dir(indexes[k],ghost):
                                                                        'not_started'})
                            # Calculate molecular weight
                            totmass = sum(cluster.mass(i) for i in range(cluster.natom()) if cluster.Z(i))
                            db[method][n]['MW'].update({n_body.ghost_dir(indexes[k],ghost):totmass})
                            # Reactivate fragments for next round
                            cluster.set_active_fragments(list(ghost))
                            cluster.update_geometry()

                    # SSFC or no BSSE
                    else:
                        cluster.set_name('{}_{}'.format(molecule.name(),n_body.cluster_dir(indexes[k])))
                        #molecule.set_name('{}_{}'.format(molecule.name(),cluster_dir(indexes[k])))
                        directory = '{}/{}/{}'.format(method,n_body.n_body_dir(n),n_body.cluster_dir(indexes[k]))
                        n_body.plant(cluster, db, kwargs, method, directory)
    
                        # Update database job_status dict
                        db[method][n]['job_status'].update({n_body.cluster_dir(indexes[k]):
                                                                       'not_started'})
                        # Calculate molecular weight
                        totmass = sum(cluster.mass(i) for i in range(cluster.natom()) if cluster.Z(i))
                        db[method][n]['MW'].update({n_body.cluster_dir(indexes[k]):totmass})
        # Check for zero jobs (happens occasionally due to cutoff)
        for method in db['methods']:
            for n in range(1, db[method]['n_body_max']+1):
                if db[method][n]['total_num_jobs'] == 0:
                    db[method]['n_body_max'] = n-1
                    break
        # Inputs successfully generated
        db['inputs_generated'] = True

    # For now open and close database frequently, until I figure out atexit
    db.close()
    db = shelve.open('database',writeback=True)

    # Check status of jobs
    if not db['jobs_complete']:
        n_incomplete = 0
        # Check all methods
        for method in db['methods'].keys():
            #### NOTE: dft_methods is defunct #####
#            if method in n_body.dft_methods:
#                outname = 'input.log'
#                complete_message = 'Normal termination of Gaussian'
#            if (method == 'b3lyp'):
            if method in n_body.dft_methods:
                outname = 'input.log'
                complete_message = 'Normal termination of Gaussian'
                error_message = 'Error termination'
            else:
                outname = 'output.dat'
                complete_message = 'Psi4 exiting successfully'
                error_message = 'Psi4 encountered an error'
            # Check all n_body_levels
            print('Before job checking:')
            for field in db[method]['farm']:
                num_fin = db[method][field]['num_jobs_complete']
                tot_num = db[method][field]['total_num_jobs']
                print('{}/{} {}-body jobs finished'.format(num_fin,tot_num,field))
                n_complete = num_fin
                if num_fin != tot_num:
                    db_stat = db[method][field]['job_status']
                    n_complete = 0
                    for job,status in db_stat.items():
                        if status == 'complete':
                            n_complete += 1
                        elif status in ('not_started', 'running', 'error'):
                            try:
                                outfile = open('{}/{}/{}/{}'.format(method,n_body.n_body_dir(field),job,outname))
                                for line in outfile:
                                    if complete_message in line:
                                        # This actually marks the job complete as soon as
                                        # the first route section is complete in the case
                                        # of g09, not when the job is totally finished.
                                        db_stat[job] = 'complete'
                                        n_complete += 1
                                        break
                                    elif error_message in line:
                                        db_stat[job] = 'error'
                                        break
                                else:
                                    db_stat[job] = 'running'
                                    n_incomplete += 1
                            except:
                                # this print statement is super annoying, muting it for now
#                                print('Exception: probably could not open file for {}'.format(job))
                                n_incomplete += 1
                db[method][field]['num_jobs_complete'] = n_complete
        if n_incomplete == 0:
            db['jobs_complete'] = True

    db.close()
    db = shelve.open('database',writeback=True)

    # Gather results
    if not db['results_computed']:
        for method in db['methods'].keys():
            print('\nAfter job checking:') 
            for field in db[method]['farm']:
                num_fin = db[method][field]['num_jobs_complete']
                tot_num = db[method][field]['total_num_jobs']
                print('{}/{} {}-body jobs finished'.format(num_fin,tot_num,field))
                if (db[method][field]['num_jobs_complete'] == db[method][field]['total_num_jobs']):
#                    if method in n_body.dft_methods:
#                        n_body.harvest_g09(db,method,field)
#                    if (method == 'b3lyp'):
                    if method in n_body.dft_methods:
                        n_body.harvest_g09(db,method,field)
                    else:
                        n_body.harvest_data(db,method,field)
            for field in db[method]['farm']:
                if (db[method][field]['num_jobs_complete'] == db[method][field]['total_num_jobs']):
#                if isinstance(field, int):
                    n_body.cook_data(db,method,field)
#                    #if db['bsse'] == 'vmfc' and field > 1:
#                    #    n_body.vmfc_cook(db, method, field)
#                    #    n_body.mbcp_cook(db, method, field)
#                    if db['bsse'] == 'vmfc':
#                        n_body.vmfc_cook(db, method, field)
#                        if field > 1:
#                            n_body.mbcp_cook(db, method, field)
#                    elif db['bsse'] == 'mbcp' and field > 1:
#                        n_body.mbcp_cook(db,method,field)
#                    # Future location of an if check for mass adjust or not
#                    if 'rotation' in db[method]['results']:
#                        n_body.mass_adjust_rotations(db,method,field)

    db.close()
    db = shelve.open('database',writeback=True)

    # Print banner to output file including user options
    n_body.banner(db)

    db.close()

    pass


##### END TAYLOR'S run_n_body() FUNCTION #####



# Integration with driver routines
#psi4.driver.procedures['energy']['n_body'] = run_n_body_energy
#psi4.driver.procedures['properties']['n_body'] = run_n_body_properties

