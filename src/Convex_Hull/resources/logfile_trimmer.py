# this script deletes logfiles not needed for CCK* from K* runs

import glob
import os

deleted_files = 0

for logdir in glob.glob("match*-kstar/kstar-*/ensembles/*csv"):

    print("Deleting %s" % logdir)

    os.remove(logdir)

    deleted_files += 1

print("Deleted %s files" % deleted_files)
