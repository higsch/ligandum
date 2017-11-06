#!/usr/bin/env python3
# encoding: utf-8

import pymzml
import os
import pickle
from ratios import Ratios
from rpy2.robjects.packages import importr

def main():
    quant_summary = '/Users/MS/Desktop/special_projects/SMHacker/quant_summary.csv'
    rt_info_file = '/Users/MS/Desktop/special_projects/SMHacker/ligand_quant_res.csv'
    
    results_class = pickle.load(
        open(
            '/Users/MS/Desktop/special_projects/SMHacker/pyQms_results.pkl',
            'rb'
        )
    )
    
    rs = Ratios(quant_summary, rt_info_file, results_class, ['TEV_H', 'TEV_L'])
    rs.read_and_parse_files()
    gen = rs.calculate_ratios('TEV_H', 'TEV_L', 'auc in window')
    for key, ratio in gen:
        print(key[0][:int(key[2])]+'*'+key[0][int(key[2]):], key[1], ratio)
    return

if __name__ == '__main__':
    main()