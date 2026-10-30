#!/bin/env python
from __future__ import absolute_import
import sys
from .lib.Database import Database
def main():
	organ, taxon = sys.argv[1:3]
	db = Database(organ=organ, taxon=taxon)
	db.checkdb(untar=False)

if __name__ == '__main__':
	main()
