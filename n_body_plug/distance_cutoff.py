import shelve
import collections
import copy

def main():
    '''
    Main function
    '''
    method = 'cam-b3lyp'
    kwargs = ['cutoff_rotation', 'cutoff_solute_rotation']
    db = shelve.open(db_name, writeback=True)
    
    if sanity_check(db):
        print("Sanity check passed!")

    extend(db, kwargs)

    n_max = db[method]['n_body_max']
    for result in kwarg:
        for n in range(1,n_max+1):
            db[method][n][kwarg]['raw_data'] = cutoff(db, method, n, 'rotation')

    # Need to generate close_solvent here I think, rather than at db generation
    # Final step: Cook cutoff data


def cutoff(db, method, n, full_result, dist = 0)
    '''
    Removes raw_data from the full_result based on the distance cutoff given and returns it
    NOTE: Need to employ a pass-able distance cutoff rather than the hardcoded database one
    '''
    cut_result = collections.OrderedDict()
    cut_result[n] = copy.deepcopy(db[method][n][full_result]['raw_data'])
    for job in list(cut_result[n].keys()):
        job_set = set(job.split('_'))
        if job_set <= set(db['close_solvent']):
            continue
        else:
            del cut_result[n][job]
    return cut_result



def extend(db, method, kwargs)
    '''
    Extends the database to hold results from kwargs list
    '''
    farm = copy.deepcopy(db[method]['farm'])
    for field in farm:
        for result in kwarg:
            db[method][field][result] = collections.OrderedDict()
            db[method][field][result]['correction'] = 0
            db[method][field][result][result] = 0
            db[method][field][result]['raw_data'] = 0
            db[method][field][result]['cooked_data'] = 0


def sanity_check(db):
    '''
    Checks a database for: 
    jobs_complete
    results_computed
    distance
    '''
    if not db['jobs_complete']:
        raise Exception("Jobs are not complete!")
    elif not db['results_computed']:
        raise Exception("Basic results have not been computed!")
    elif not db['distance']:
        raise Exception("No distance cutoff given!")
    else:
        return True


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        db_name = str(sys.argv[1])
    else: 
        db_name = 'database'
    main()
