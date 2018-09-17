import shelve
import collections
import copy

def main():
    '''
    Main function
    '''
    # Define method, result(s) to add, grab database 
    method = 'cam-b3lyp'
    results = ['cutoff_rotation', 'cutoff_solute_rotation']
#    distance = 5
    distance = dist_kwrg
    db = shelve.open(db_name, writeback=True)
    n_max = db[method]['n_body_max']
    
    # Make sure database has been filled out by n_body driver
    if sanity_check(db):
        print("Sanity check passed!")
    else: # Should raise Exception for issue
        print("Sanity check not passed!")

    # NOTE: Need to generate close_solvent here I think, rather than at db generation
    for result in results:
        # Extend the database to hold the new result(s)
        db = extend(db, method, result, distance)

        # Cut the irrelevant data for result(s) out of 'raw_data' and store
        for n in range(1,n_max+1):
            db[method][n][result]['raw_data'] = cutoff(db, method, n, 
                                                       result.replace('cutoff_',''))

            # Cook new result(s) data
            cook(db, method, n, result)


def cook(db, method, n, result):
    '''
    Cook the data as in the n_body driver, one result at a time
    '''
    # Automatic cooking requires all result data be stored as lists in the
    # database. Even if the result is a single number (i.e. energies)
    cooked_data = collections.OrderedDict()
    raw_data    = collections.OrderedDict()

    # Read in all cooked_data for m < n
    for m in range(1, n):
        cooked_data[m] = copy.deepcopy(db[method][m][result]['cooked_data'])

    # Read in raw data for m
    cooked_data[n] = copy.deepcopy(db[method][n][result]['raw_data'])
    if len(cooked_data[n]) == 0: # No raw_data, cutoff probably cut all n-body jobs
        print('No raw_data for {}-body {}, distance_cutoff has likely been reached'.format(n,result))
        db[method][n][result][result] = copy.deepcopy(db[method][n-1][result][result])
        return

    # Cook n data
    for m in range(1,n):
        for n_key in cooked_data[n]:
            for m_key in cooked_data[m]:
                m_set = set(m_key.split('_'))
                n_set  = set(n_key.split('_'))
                if len(m_set.intersection(n_set)) == len(m_set):
                    cooked_data[n][n_key] = [x-y for x,y in 
                    zip(cooked_data[n][n_key],cooked_data[m][m_key])]

    # Set correction list length
    correction = [0.0 for i in range(len(next(iter(cooked_data[n].values()))))]

    # Add up all contributions to current correction
    for key,val in cooked_data[n].items():
        correction = [x+y for x,y in zip(correction,val)]

    db[method][n][result]['cooked_data'] = cooked_data[n]
    db[method][n][result]['correction']  = copy.deepcopy(correction)

    # Final result value is the correction for n plus the result value from n-1
    db[method][n][result][result] = copy.deepcopy(correction)
    if n > 1: 
        prev = db[method][n-1][result][result]
        curr = db[method][n][result][result]
        if prev != 0:
            db[method][n][result][result] = [x+y for x,y in zip(curr,prev)]
        else:
            print("Previous result 0, {}-body jobs are likely not finished. Stopping harvest.".format(n-1))

    # Wipe out data dicts for next result
    cooked_data.clear()
    raw_data.clear()


def cutoff(db, method, n, full_result):
    '''
    Removes raw_data from the full_result based on the distance cutoff given and returns it
    NOTE: Need to employ a pass-able distance cutoff rather than the hardcoded database one
    '''
    cut_result = copy.deepcopy(db[method][n][full_result]['raw_data']) # Copy full raw_data
    for job in list(cut_result.keys()):
        job_set = set(job.split('_'))
        if job_set <= set(db['close_solvent']): # Keep if all frags in close_solvent
            continue
        else:
            del cut_result[job] # Dump if any frags outside close_solvent
    return cut_result



def extend(db, method, result, dist = False):
    '''
    Extends the database to hold result
    Adds distance cutoff to db if not already present
    Determines close_solvent molecules
    '''
    # Update distance
    db['distance'] = dist

    # Determine close solvent molecules
    if dist:
        db['close_solvent'] = []
        for slt in db['solvent_distances']:
            if db['solvent_distances'][slt] < db['distance']:
                db['close_solvent'].append(slt)

    # Extend database to hold result
    if not result in db[method][1]: # Dict should be properly extended if 1 exists
        farm = copy.deepcopy(db[method]['farm']) # All n-body
        for field in farm:
            db[method][field][result] = collections.OrderedDict()
            db[method][field][result]['correction'] = 0
            db[method][field][result][result] = 0
            db[method][field][result]['raw_data'] = collections.OrderedDict()
            db[method][field][result]['cooked_data'] = collections.OrderedDict()

    return db


def sanity_check(db):
    '''
    Checks a database for: 
    jobs_complete
    results_computed
    close_solvent
    '''
    if not db['jobs_complete']:
        raise Exception("Jobs are not complete!")
    elif not db['results_computed']:
        raise Exception("Basic results have not been computed!")
    elif not db['solvent_distances']:
        raise Exception("Close solvent molecules not determined by the driver!")
    else:
        return True


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        db_name = str(sys.argv[1])
    else: 
        db_name = 'database'
    if len(sys.argv) == 3:
        dist_kwrg = float(sys.argv[2])
    main()
