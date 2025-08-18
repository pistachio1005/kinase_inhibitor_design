# this helper script checks search progress for all search in pair_set_kstar

import glob

total_matches = 0
run_match_count = 0
done_match_count = 0
completed_matches = []
running_matches = []

for matchdir in glob.glob("FSV_Scan2/match*-kstar"):

    total_matches += 1

    print("\n--Checking %s--" % matchdir)

    total_doublets = 0
    running_doublets = []
    finished_doublets = 0

    doublet_dir = ("%s/%s" % (matchdir, "kstar-*"))

    for d in glob.glob(doublet_dir):
        total_doublets += 1

    logfile = ("%s/%s" % (doublet_dir, "submit.out"))

    for log in glob.glob(logfile):
        is_done = False
        for line in open(log, 'r'):
            if 'completed' in line:
                is_done = True
        if is_done:
            finished_doublets += 1
        if not is_done:
            running_doublets.append(log.split('/')[2])

    if finished_doublets > total_doublets:
        print("ERROR! Something went wrong with job checking")

    elif finished_doublets == total_doublets:
        print("This match has finished all %s doublets!" % total_doublets)
        done_match_count += 1
        completed_matches.append(matchdir.split('/')[1])

    elif finished_doublets < total_doublets:
        print("%s of %s doublets have finished running" % (finished_doublets, total_doublets))
        print("Still running: %s" % running_doublets)
        run_match_count += 1
        running_matches.append(matchdir.split('/')[1])

print("\n---- Overall stats ----\n")
print("%s total matches\n" % total_matches)
print("%s matches still running:" % run_match_count)
print(running_matches)
print("\n%s matches finished:" % done_match_count)
print(completed_matches)
