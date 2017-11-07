#!/usr/bin/env python3
# encoding: utf-8

import sys
import bisect
import pyqms
import codecs
import csv
import math
from collections import namedtuple
from rpy2.robjects.lib.grdevices import grdevices
from numpy import empty

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
            result = 0.0
            
            try:
                float1 = float(value[label1][quant_field])
            except (KeyError, ValueError, TypeError) as e:
                #print(e)
                float1 = 0.0
            
            try:
                float2 = float(value[label2][quant_field])
            except (KeyError, ValueError, TypeError) as e:
                #print(e)
                float2 = 0.0
                
            try:
                result = float1/float2
            except ZeroDivisionError as e:
                #print(e)
                result = 20.0
            except:
                print('Unknown error in ratio calculation!')
            
            #print(float1, float2)
            if float1 == 0.0 and float2 == 0.0:
                yield key, 0.0
                continue
            
            yield key, result
    
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
    
    
    def _import_rpy2(self):
        
        try:
            global rpy2
            global importr
            global robjects
            global r
            global graphics
            global grdevices
            import rpy2
            from rpy2.robjects.packages import importr
            import rpy2.robjects as robjects
            from rpy2.robjects import r
            graphics = importr('graphics')
            grdevices = importr('grDevices')
            success = True
        except:
            success = False
        return success
    
    
    def _init_r_plot(self, file_name):
        assert self._import_rpy2() is True, 'require R & rpy2 installed...'
        graphics = importr('graphics')
        grdevices = importr('grDevices')
        grdevices.pdf(file_name)
        return graphics, grdevices
    
    
    def _generate_r_colors(self, tag, elements):
        p = {
            'rainbow' : list(r.rainbow( elements , start=0.2, end=1 )),
            'terrain' : list(r('{0}.colors'.format( 'terrain' ))( elements )),
            'topo'    : list(r('{0}.colors'.format( 'cm' ))( elements )),
            'heat'    : list(r('{0}.colors'.format( 'heat' ))( elements )),
        }
        reverse = False
        if tag[-2:] == '_r':
            reverse = True
            tag = tag[:-2]
        if tag not in p.keys():
            print('Do not {0} as color tag, falling back to rainbow')
            tag = 'rainbow'
        colors = p[ tag ]
        if reverse:
            colors.reverse()
        return colors
    
    
    def colorize_score(self, score ):
        
        color = [0, 0, 0]
        colorGradient = [ (score_threshold, rgb_tuple) for (score_threshold, rgb_tuple) in sorted(pyqms.params['COLORS'].items()) ]
        if score is not None:
            idx = bisect.bisect(  colorGradient, ( score, )  )
            if idx == 0:
                color =  colorGradient[0][1]
            elif idx == len( colorGradient):
                color =  colorGradient[-1][1]
            else:
                # linear interpolation ... between idx-1 & idx
                dX = ( score -  colorGradient[ idx - 1 ][ 0 ] ) / (  colorGradient[ idx ][ 0 ] - colorGradient[ idx - 1 ][ 0 ] )
                for color_chanel in range(3):
                    d_ = dX * ( colorGradient[ idx ][ 1 ][ color_chanel ] -  colorGradient[ idx - 1 ][ 1 ][ color_chanel ])
                    if abs( d_ ) <= sys.float_info.epsilon :
                        color[ color_chanel ] = int(round( colorGradient[idx - 1][ 1 ][ color_chanel ]))
                    else:
                        color[ color_chanel ] = int(round( colorGradient[idx - 1][ 1 ][ color_chanel ] + d_))
        hexed_color = '#' + ''.join( [ hex(c)[2:] for c in color ] )

        return color, hexed_color


    def plot_pairs(
            self,
            key_list,
            file_name,
            rt_offset          = None,
            xlimits            = None,
            title              = None,
            grdevices          = None,
            graphics           = None
                ):

        zlimits = [0, 1]
        
        if rt_offset is None:
            rt_offset = 3.0
        
        colors = []
        
        for n in range(0, 101, 1):
            rgb_col, hex_col =  self.colorize_score(n/100)
            colors.append(hex_col)
        
        for n, key in enumerate(key_list):
            
            if key not in self.keys():
                print('Warning, do not have match results for {0}'.format(key))
                continue
            
            ms2_evidences = {}
            
            curve_dict = {}
            for label_key, label_value in self[key].items():
                
                if label_value['has_MS2_id'] == True:
                    ms2_evidences.update({label_key: label_value['evidences (min)']})
                
                x = []
                y = []
                s = []
                c = []
                for entry in label_value['data']:
                    if float(label_value['start (min)']) - rt_offset <= entry.rt <= float(label_value['stop (min)']) + rt_offset:
                        x.append(entry.rt)
                        y.append(entry.scaling_factor)
                        s.append(entry.score)
                        
                if x == []:
                    continue
                
                for score in s:
                    c.append(colors[int(round(score*100))])
            
                if xlimits is None:
                    if len(x) >= 1:
                        xlimits = [x[0], x[-1]]
                    
                curve_dict.update({
                        label_key: {
                                'x': x,
                                'y': y,
                                's': s,
                                'c': c,
                                'max_y': max(y),
                                'xlimits': xlimits
                            }
                    })
            
            # determine overall max_y
            max_list = []
            min_xlimits_list = []
            max_xlimits_list = []
            for label_key, label_value in curve_dict.items():
                max_list.append(label_value['max_y'])
                if label_value['xlimits'] is not None:
                    min_xlimits_list.append(label_value['xlimits'][0])
                    max_xlimits_list.append(label_value['xlimits'][-1])
            
            if max_list == [] or min_xlimits_list == []:
                continue
            
            if grdevices is None or graphics is None:
                graphics, grdevices = self._init_r_plot(file_name)
            
            graphics.par(mfrow = robjects.FloatVector([len(key_list), 1]),)
            
            max_y = max(max_list)
            min_x_limit = min(min_xlimits_list) - 1
            max_x_limit = max(max_xlimits_list) + 1
            
            params = {
                'pch' : 19,
                'cex' : 0.7,
                'xlab' : 'Retention Time [min]',
                'ylab' : 'Relative  Abundance',
                'xlim' : r.c(min_x_limit, max_x_limit),
                'ylim' : r.c(0, max_y * 1.1),
                'main' : '{0}\n Charge: {1}\n Mods: {2}'.format(key[0][:int(key[2])]+'*'+key[0][int(key[2]):], key[1], key[3]),
                'frame': False
            }
            
            graphics.plot(
                robjects.FloatVector([]),
                robjects.FloatVector([]),
                **params
            )
            
            k = 0
            for key, value in curve_dict.items():
                graphics.lines(
                    robjects.FloatVector(value['x']),
                    robjects.FloatVector(value['y']),
                    type = 'l',
                    lwd = 0.7,
                    col = grdevices.palette()[k]
                )
                    
                graphics.points(
                    robjects.FloatVector(value['x']),
                    robjects.FloatVector(value['y']),
                    col = robjects.StrVector(value['c']),
                    lwd = 0.1,
                    pch = 19
                )
                
                try:
                    evidences = self._parse_evidences(ms2_evidences[key])
                    if len(evidences) > 0:
                        graphics.points(
                            robjects.FloatVector(evidences),
                            robjects.FloatVector([0]*len(evidences)),
                            col = grdevices.palette()[k],
                            lwd = 0.1,
                            pch = 24
                        )
                except KeyError:
                    pass
                
                k += 1
                
            
            
            
            grdevices.dev_off()

        return grdevices
    
    
    def _parse_evidences(self,
            evidences
                ):
        
        evidences_list = []
        for evidence in evidences.split(';'):
            evidences_list.append(evidence[evidence.find('@')+1:])
            
        return evidences_list