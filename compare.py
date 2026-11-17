from __future__ import absolute_import
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
sys.path.insert(0, str(here))

from lib.Comparative import main as _main

def main():
    _main()

if __name__ == '__main__':
    main()

