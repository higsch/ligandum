#!/usr/bin/env python3
# encoding: utf-8

import pickle
import os
from _operator import xor

def main():
    pickle_file = '/Users/MS/Desktop/special_projects/SMHacker/pyQms_results.pkl'
    results = pickle.load(
        open( pickle_file, 'rb')
    )
    out_folder = os.path.join(
        os.path.dirname(pickle_file),
        'plots'
    )
    if os.path.exists(out_folder) is False:
        os.mkdir(out_folder)
    print('Plotting into folder: {0}'.format(out_folder))
    for n, key in enumerate(results.keys()):
        if key[1] != 'C(112)H(179)N(33)O(31)' and key[1] != 'C(107)H(179)13C(5)15N(1)N(32)O(31)':
            continue
#         if len(results[key]['data']) <= 10:
#             continue
        mzml_filename = key.file_name
        if os.sep in mzml_filename:
            mzml_filename = os.path.basename(mzml_filename)

        file_name = os.path.join(
            out_folder ,
            'MIC_2D_{0}_{1}_{2}_{3}.pdf'.format(
                '_'.join(
                    results.lookup['formula to molecule'][ key.formula ]
                ),
                key.charge,
                key.label_percentiles,
                mzml_filename
            )
        )
        graphics, grdevices = results.init_r_plot(file_name)
        results.plot_MICs_2D(
            [key],
            graphics = graphics
        )


    return


if __name__ == '__main__':
    main()
