#!/usr/bin/env python3
# encoding: utf-8

import pymzml
import os
from rpy2.robjects.packages import importr

def main():
    print(os.environ.get("R_HOME"))
    run = pymzml.run.Reader('/Users/MS/Desktop/special_projects/SMHacker/170209_SMH_170205_P9_05_new.mzML')

    for spectrum in run:
        if spectrum['ms level'] == 1:
            print(spectrum)
            break
    return

if __name__ == '__main__':
    main()