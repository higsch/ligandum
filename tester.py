#!/usr/bin/env python3
# encoding: utf-8

import pickle
from ratios import Ratios

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
#     for key, value in rs.get_results_by_sequence('EVDEQMLNVQNKNSSYFVEWIPNNVK'):
#         print(key, value)
#     gen = rs.calculate_ratios('TEV_H', 'TEV_L', 'auc in window')
#     for key, ratio in sorted(gen):
#         print(key[0][:int(key[2])]+'*'+key[0][int(key[2]):], key[1], ratio)
    
    rs.curate_pairs()
    
    for key, value in rs.get_results_by_sequence('ATFSEFAAKHAK', labels_only = False):
         print(key, value)
    
#     gen = rs.calculate_ratios('TEV_H', 'TEV_L', 'max I in window')
#     i = 0
#     for key, ratio in sorted(gen):
#         rs.plot_pairs([key], '/Users/MS/Desktop/special_projects/SMHacker/plots/plot_{0}_{1}.pdf'.format(key[0], key[1]), {'TEV_H': 0, 'TEV_L': 1})
#         i += 1
#         if i % 20 == 0:
#             print('> Plot tuples:', i, end = '\r')
    

    return

if __name__ == '__main__':
    main()