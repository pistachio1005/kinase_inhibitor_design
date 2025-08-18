# for extracting GMEC PDBs after ABAA

import os
import glob
import shutil

out_directory = "ABAA_GMEC"

os.mkdir(out_directory)

print("Storing outputs in %s" % out_directory)

# check the logfile and find out which runs were successful (> 0 K* score, positive partition functions)
bad_runs = 0
good_runs = []
total_runs = 0
for logfile in glob.glob("ABAA/match*/submit.out"):

    total_runs += 1

    matchname = logfile.split('/')[1].split('-')[0]
    onheader = False

    f = open(logfile, 'r')
    for line in f:

        if onheader:

            log_fields = line.split(',')[2:]
            is_bad_run = False

            for score in log_fields:
                if 'none' in score or '-' in score:
                    print("SKIP: Negative or 0 K* score or pfunc for %s!" % logfile)
                    bad_runs += 1
                    is_bad_run = True
                    onheader = False
                    break

            if not is_bad_run:
                print("Positive K* score and pfuncs for %s" % logfile)
                good_runs.append(matchname)
                onheader = False

        if 'Assignments' in line:
            onheader = True

print("\n%s out of %s matches had 0 or negative K* score. These will be ignored.\n" % (bad_runs, total_runs))

# copy over all the PDB files from the K* runs (in ABBA/)
for pdb in glob.glob("ABAA/match*/ensembles/*pdb"):

    match_name = pdb.split('/')[1].split('-')[0]

    if match_name in good_runs:
        new_pdb_name = pdb.split('/')[1].split('-')[0] + "-ABAA.pdb"
        new_path = os.path.join("ABAA_GMEC", new_pdb_name)
        shutil.copy(pdb, new_path)
        print("copied file %s" % new_pdb_name)
    else:
        print("Skipping match %s" % match_name)

print("Copied over %s good runs to %s" % (total_runs - bad_runs, out_directory))

# print out how ABAA rank by K* score (can be helpful for prioritizing which matches to run)
all_matches_tracker = {}

for logfile in glob.glob("ABAA/match*/submit.out"):

    matchname = logfile.split('/')[1]

    f = open(logfile, 'r')
    on_start = False
    for line in f:
        if on_start:
            score = line.split(',')[2]
            if score == 'none':
                score = 0
            else:
                score = float(score)
            all_matches_tracker[matchname] = score
            break
        if 'Seq' in line:
            on_start = True

print("\nABAA results ordered by K* score:")
ordered_match_tracker = {k: v for k, v in reversed(sorted(all_matches_tracker.items(), key=lambda item: item[1]))}
for k, v in ordered_match_tracker.items():
    print(k, v)
