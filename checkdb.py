import sys
from lib.Database import Database
organ, taxon = sys.argv[1:3]
db = Database(organ=organ, taxon=taxon)
db.checkdb(untar=False)

