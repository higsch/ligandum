#!/usr/bin/env python3
# encoding: utf-8

import pymzml
import os
from ratios import Ratios
from rpy2.robjects.packages import importr

def main():
    quant_summary = '/Users/MS/Desktop/special_projects/SMHacker/quant_summary.csv'
    rt_info_file = '/Users/MS/Desktop/special_projects/SMHacker/ligand_quant_res.csv'
    
    rs = Ratios(quant_summary, rt_info_file, ['TEV_H', 'TEV_L'])
    rs.read_and_parse_files()
    print(rs.keys())
    return

if __name__ == '__main__':
    main()