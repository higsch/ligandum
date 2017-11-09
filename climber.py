#!/usr/bin/env python3
# encoding: utf-8

default_params = {
        'INTENSITY_PERCENTILE': 0.5,
        'RT_OFFSET_TO_MS2'    : 1
    }

class Climber:
    
    def __init__(
            self,
            params = None
                ):
        
        if params is None:
            self.params = default_params
        else:
            self.params = params
        
        pass
    return