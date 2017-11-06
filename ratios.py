#!/usr/bin/env python3
# encoding: utf-8

from __future__ import division
import pyqms
import codecs
import csv
import numbers
from collections import namedtuple

r_key = namedtuple(
    'r_key',
    [
        'sequence',
        'charge',
        'label_position',
        'mods'
    ]
)

match = namedtuple(
    'match',
    [
        'spec_id',
        'rt',
        'score',
        'scaling_factor',
        'peaks'
    ]
)

default_body_fields = [
        'file_name',
        'trivial_name(s)',
        'evidences (min)',
        'formula',
        'molecule',
        'max I in window (rt)',
        'stop (min)',
        'start (min)',
        'auc in window',
        'auc in window',
        'max I in window',
        'sum I in window',
        'max I in window (score)',
        'label_percentiles'
    ]


class Ratios(dict):
    
    def __init__(
            self,
            quant_summary_file = None,
            rt_info_file       = None,
            results_class      = None,
            labels             = None
                ):
        
        self.quant_summary_file = quant_summary_file
        self.rt_info_file = rt_info_file
        self.results_class = results_class
        self.labels = labels
        
        self._r_key_class = r_key
        self._match_class = match
        
        self.quant_summary_dicts = None
        self.rt_info_dicts = None
        
        return
    
    
    def calculate_ratios(
                self,
                label1,
                label2,
                quant_field
                    ):
        
        for key, value in self.items():
            
            if label1 not in value.keys():
                continue
            if label2 not in value.keys():
                continue
            try:
                yield key, float(value[label1][quant_field])/float(value[label2][quant_field])
            except TypeError:
                print('no number!')
            except ZeroDivisionError:
                yield key, 'zero division'
            except ValueError:
                yield key, 'n.d.'
    
    
    def add_body(
            self, 
            key, 
            label,
            info
                ):
        
        r_key = self._r_key_class(*key)
        
        if r_key in self:
            if label not in self[r_key]:
                self[r_key].update({label: info})
        else:
            self[r_key] = {
                    label: info
                }
        
        return r_key
    
    
    def read_and_parse_files(
            self,
            quant_summary_file = None,
            rt_info_file       = None,
            results_class      = None,
            labels             = None
                ):
        
        print('Ratios are prepared...')
        
        if quant_summary_file is not None:
            self.quant_summary_file = quant_summary_file
        if rt_info_file is not None:
            self.rt_info_file = rt_info_file
        if results_class is not None:
            self.results_class = results_class
        if labels is not None:
            self.labels = labels
        
        
        if self.rt_info_file is not None:
            self.rt_info_dicts = self._get_dicts_from_file(self.rt_info_file)
             
        if self.quant_summary_file is not None:
            self.quant_summary_dicts = self._get_dicts_from_file(self.quant_summary_file)
        
        # read summary file dict in common data structure
        for line in self.quant_summary_dicts:
            molecule, label, label_position, mods = self._split_molecule(line['molecule'])
            self.add_body(
                    key = (molecule, line['charge'], label_position, mods),
                    label = label,
                    info = self._extract_molecule_info(line)
                )
        
        # match results class to quant summary
        print('Evaluating {0} peptide tuples...'.format(len(self)))
        for result_key, result_value in self.items():
            #print(result_key[0])
            for label_key, label_value in result_value.items():
                for key, i, entry in self.results_class.extract_results(
                        molecules         = None,
                        charges           = [int(result_key[1])],
                        file_names        = label_value['file_name'],
                        label_percentiles = None,
                        formulas          = label_value['formula'],
                        score_threshold   = None
                    ):
                    self[result_key][label_key]['data'].append(entry)
        return


    def _get_dicts_from_file(
            self, 
            file
                ):
        
        tmp_dicts = []
        if file.endswith('.csv'):
            with codecs.open(file, mode='r', encoding='utf-8'  ) as rif:
                dict_reader = csv.DictReader(rif)
                for line_dict in dict_reader:
                    tmp_dicts.append(line_dict)
        elif file.endswith('.xlsx'):
            tmp_dicts = pyqms.adaptors.read_xlsx_file(file)
        else:
            print('Extension: {0} of file {1} not recognized'.format(file.split('.')[-1], file))
            exit(1)
        
        return tmp_dicts
    
    
    def _split_molecule(
            self,
            molecule
                ):
    
        label = None
        label_position = -1
        mod_start = molecule.find('#')
        for label_name in self.labels:
            pos = molecule.find(label_name)
            if pos >= 0:
                label = label_name
                sub_string_label = molecule[pos+len(label_name)+1:]
                pos_end_label = sub_string_label.find(';')
                if pos_end_label >= 0:
                    label_position = sub_string_label[:sub_string_label.find(';')]
                else:
                    label_position = sub_string_label
                break
        
        mods = molecule[mod_start+1:].replace('{0}:{1}'.format(label, label_position), '')
        mods = mods.replace(';;', ';')
        if mods[0:1] == ';':
            mods = mods[1:]
        elif mods[-1:] == ';':
            mods = mods[:-1]
            
        return molecule[:mod_start], label, label_position, mods
    
    
    def _extract_molecule_info(
            self, 
            line
        ):
        
        tmp = {}
        for field_name in default_body_fields:
            tmp[field_name] = line[field_name]
        
        tmp['data'] = []
        
        if tmp['trivial_name(s)'].find('no MS2;') >= 0:
            tmp['has_MS2_id'] = False
            tmp['trivial_name(s)'] = tmp['trivial_name(s)'].replace('no MS2;', '')
        else:
            tmp['has_MS2_id'] = True
        
        return tmp
    