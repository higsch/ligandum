#!/usr/bin/env python3
# encoding: utf-8
'''
Ligandum -- Analysis of ligandability mass spectrometry

Created by Matthias Stahl at TU Muenchen, 2017
'''

import ursgal
import pymzml
import pyqms
import sys
import pickle
import os
import pprint
import random


def showStartHello():
    print('{0:-^100}'.format('###'))
    print('{0:-^100}'.format('This is Ligandum speaking! Let\'s start...'))
    print('{0:-^100}'.format('###'))
    print()
    return


def new_userdefined_unimod_molecule(mass, name, composition):
    # make modification dictionary
    modification_dict = {
        'mass' : mass,
        'name' : name,
        'composition' : composition,
        'id': '1'
    }
    
    # initialize unimod mapper and write to unimod.xml
    u = ursgal.UnimodMapper()
    u.writeXML(modification_dict, xmlFile = '/Library/Frameworks/Python.framework/Versions/3.5/lib/python3.5/site-packages/pyqms/kb/ext/unimod.xml')
    return

    
def msms_identification(mzml_file, database_file):
    # define parameters
    params = {
        'enzyme': 'trypsin',
        'frag_mass_tolerance': 0.5,
        'frag_mass_tolerance_unit': 'da',
        'decoy_generation_mode' : 'reverse_protein',
        'precursor_min_charge' : '2',
        'modifications' : [
            'M,opt,any,Oxidation',
            'C,fix,any,Carbamidomethyl',
            'K,opt,any,TEV_H',
            'K,opt,any,TEV_L'
        ],
    }
    
    # initialize Ursgal controller
    uc = ursgal.UController(
        params = params
    )
    
    # generate reverse protein sequences and initialize new database
    new_target_decoy_db_name = uc.generate_target_decoy(
        input_files = database_file,
        output_file_name = 'new_target_decoy.fasta',
    )
    print('Generated target decoy database: {0}'.format(new_target_decoy_db_name))
    uc.params['database'] = new_target_decoy_db_name
    
    # perform search with search engine (writes output files to file system)
    search_result = uc.search(
        input_file = mzml_file,
        engine = 'msgfplus_v2016_09_16'
    )
    
    # validate search engine results with percolator (writes output files to file system)
    validated_result = uc.validate(
        input_file = search_result,
        engine     = 'percolator_2_08',
    )
    
    return validated_result


def ligandability_quantification(mzml_file, molecule_list, evidence_lookup, formatted_fixed_labels):
    run = pymzml.run.Reader(mzml_file)
    params = {
        'molecules'        : molecule_list,
        'charges'          : [1, 2, 3, 4, 5],
        'fixed_labels'     : formatted_fixed_labels,
        'verbose'          : True,
        'evidences'        : evidence_lookup
    }
 
    lib = pyqms.IsotopologueLibrary( **params )
      
    results = None
    mzml_file_basename = os.path.basename(mzml_file)
    for spectrum in run:
        scan_time, unit = spectrum.get('scan time', (None, None ))
        print(unit)
        if unit == 'second':
            time_div_factor = 60.0
            # convert seconds to minutes...
        else:
            time_div_factor = 1
            
        if spectrum['ms level'] == 1:
            results = lib.match_all(
                mz_i_list = spectrum.centroidedPeaks,
                file_name = mzml_file_basename,
                spec_id   = spectrum['id'],
                spec_rt   = spectrum['scan time'] / time_div_factor,
                results   = results
            )
    return results

def edit_molecule_list(molecule_list, labels):
    # delete molecules with more than one TEV modification
    for molecule in molecule_list[:]:
        num_of_labels = 0
        for label in labels:
            num_of_labels += molecule.count(label['name'])
        if not num_of_labels == 1:
            molecule_list.remove(molecule)
    print('Initially found {0} peptides via MS2 sequencing. Now I\'ll look for nonsequenced partners.'.format(len(molecule_list)))
    
    for molecule in molecule_list[:]:
        check_pairs(molecule, molecule_list, labels)
    
    print('There are {0} peptides after partner generation.'.format(len(molecule_list)))
    return


def check_pairs(molecule, molecule_list, labels):
    partner_label_name = ''
    current_label_name = ''
    for label in labels:
        if molecule.count(label['name']) == 0:
            partner_label_name = label['name']
        else:
            current_label_name = label['name']
    
    partner_molecule = molecule.replace(current_label_name, partner_label_name)
     
    if not partner_molecule in molecule_list:
        molecule_list.append(partner_molecule)
    
    return


def main():
    showStartHello()
    
    # add new modifications
    labels = [
        {
            'name': 'TEV_H',
            'mass': '470.26338',
            'composition': {'C': 15, '13C': 5, 'H': 32, 'N': 7, '15N': 1, 'O': 5}
        },
        {
            'name': 'TEV_L',
            'mass': '464.24957',
            'composition': {'C': 20, 'H': 32, 'N': 8, 'O': 5}
        }
    ]
    
    for label in labels:
        new_userdefined_unimod_molecule(label['mass'], label['name'], label['composition'])
    
    # MS/MS identification and validation, output is written to file system
    database_file = '/Users/MS/Desktop/special_projects/SMHacker/28092017human.fasta'
    # mzml_file = '/Users/MS/Desktop/special_projects/SMHacker/170209_SMH_170205_P9_05_ultrashort.mzML'
    mzml_file = '/Users/MS/Desktop/special_projects/SMHacker/170209_SMH_170205_P9_05_short.mzML'
    validated_result = msms_identification(mzml_file, database_file)
    
    # MS isotopic ligandability quantification
    # evidence_file = '/Users/MS/Desktop/special_projects/SMHacker/msgfplus_v2016_09_16/170209_SMH_170205_P9_05_ultrashort_msgfplus_v2016_09_16_pmap_unified_percolator_validated.csv'
    evidence_file = '/Users/MS/Desktop/special_projects/SMHacker/msgfplus_v2016_09_16/170209_SMH_170205_P9_05_short_msgfplus_v2016_09_16_pmap_unified_percolator_validated.csv'
    out_folder = '/Users/MS/Desktop/special_projects/SMHacker/msgfplus_v2016_09_16'
    
    tmp_fixed_labels = {
        'C' : [
            {
                'element_composition' : {'O': 1, 'H': 3, '14N': 1, 'C': 2},
                'evidence_mod_name': 'Carbamidomethyl'
            }
        ]
    }
    formatted_fixed_labels, evidence_lookup, molecule_list = pyqms.adaptors.parse_evidence(
        fixed_labels         = tmp_fixed_labels,
        evidence_files       = [ evidence_file ],
        evidence_score_field = 'PEP'
    )
    edit_molecule_list(molecule_list, labels)
    
    results = ligandability_quantification(mzml_file, molecule_list, evidence_lookup, formatted_fixed_labels)
    results.write_result_csv(out_folder + '/ligand_quant_res.csv')
    
    # serialize, not really necessary...
    pickle.dump(
        results,
        open(
            '/Users/MS/Desktop/special_projects/SMHacker/pyQms_results.pkl',
            'wb'
        )
    )
    
    # deserialize
    results_class = pickle.load(
        open(
            '/Users/MS/Desktop/special_projects/SMHacker/pyQms_results.pkl',
            'rb'
        )
    )
    rt_border_tolerance = 5

    quant_summary_file  = '/Users/MS/Desktop/special_projects/SMHacker/quant_summary.xlsx'
    
    results_class.write_rt_info_file(
        output_file         = quant_summary_file,
        list_of_csvdicts    = None,
        trivial_name_lookup = None,
        rt_border_tolerance = rt_border_tolerance,
        update              = True
    )
    
    results_class.calc_amounts_from_rt_info_file(
        rt_info_file         = quant_summary_file,
        rt_border_tolerance  = rt_border_tolerance,
        calc_amount_function = None
    )

    return
    

if __name__ == '__main__':
    main()
