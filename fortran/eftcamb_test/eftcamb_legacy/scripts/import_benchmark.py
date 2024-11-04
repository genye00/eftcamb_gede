#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2019 by the EFTCAMB authors
#
# The EFTCAMB code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file eftcamb/LICENSE at
# the top level of the EFTCAMB distribution.
#
#----------------------------------------------------------------------------------------

"""
This python script takes care of importing and creating a dictionary with the
benchmark results.

Developed by Marco Raveri (mraveri@uchicago.edu) for the EFTCAMB code.
"""

# ***************************************************************************************

__version__ = '0.0' #: version of the application

# ***************************************************************************************

"""
import of modules
"""

import numpy as np
import math
import sys
import os

# ***************************************************************************************

"""
Hard coded options
"""

n_sigma       = 3.0       # number of sigma in time measurements to distinguish things

# ***************************************************************************************

# define the paths:
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir   = script_dir+'/../benchmark_results'
out_dir    = script_dir+'/../benchmark_plots'

# some objects to hold the benchmark results:
class EFTCAMB_benchmark:

    def __init__( self ):
        # initialize everything to zero:
        self.identifier   = 0
        self.year         = 0
        self.month        = 0
        self.day          = 0
        self.hour         = 0
        self.minute       = 0
        self.models_bench = {}
        self.models_rel_b = {}

    def is_close( self, other_bench ):
        # two benchmark results are similar if they have different identifiers
        # models_bench is within error bars.

        # check whether the two benchmarks have the same models:
        if ( sorted(self.models_bench.keys()) != sorted(other_bench.models_bench.keys()) ):
            return False
        # check whether all measurements are compatible:
        for model_name in sorted( self.models_bench.keys() ):
            value_1 = self.models_bench[model_name][0]
            err_1   = self.models_bench[model_name][1]
            value_2 = other_bench.models_bench[model_name][0]
            err_2   = other_bench.models_bench[model_name][1]
            diff_significance = abs(value_1-value_2)/math.sqrt( err_1**2 +err_2**2 )
            if diff_significance>n_sigma:
                return False
        # if no significant difference is found return true:
        return True

    def compute_relative_bench( self ):
        # get the value and error in GR:
        value_GR = self.models_bench[ '1_EFT_GR' ][0]
        err_GR   = self.models_bench[ '1_EFT_GR' ][1]
        # loop over models:
        for model in self.models_bench:
            # get value and error of the model:
            value_model = self.models_bench[ model ][0]
            err_model   = self.models_bench[ model ][1]
            # propagate errors:
            value_rel = value_model/value_GR
            err_rel   = value_rel*math.sqrt( (err_GR/value_GR)**2 +(err_model/value_model)**2 )
            # save results:
            self.models_rel_b[ model ] = [ value_rel, err_rel ]

# import the files as a nested dictionary:
benchmark_logs_files = [ file_name for file_name in os.listdir(data_dir) if 'log' in file_name ]
benchmark_results    = {}

for file_name in benchmark_logs_files:
    # split the file name:
    split_name = file_name.split(".")
    # test:
    if split_name[0]=='benchmark' and split_name[2]=='log' and len(split_name[1])==13:
        # initialize the dictionary:
        results    = EFTCAMB_benchmark()
        # get the identifier, date and time:
        results.identifier = split_name[1]
        results.year       = split_name[1][0:4]
        results.month      = split_name[1][4:6]
        results.day        = split_name[1][6:8]
        results.hour       = split_name[1][9:11]
        results.minute     = split_name[1][11:13]
        # import the data:
        with open(data_dir+'/'+file_name, 'r') as infile:
            for line in infile:
                split_line = line.split()
                results.models_bench[ split_line[0] ] = [float(split_line[1]), float(split_line[3])]
        # check the data length:
        if ( len(results.models_bench.keys())<2 ):
            print 'Rejecting file: ', file_name
            print 'Number of models < 2.'
            continue
        # check for time drifting:
        GR1 = results.models_bench['base_params']
        GR2 = results.models_bench['1_EFT_GR']
        drift_significance = abs(GR1[0]-GR2[0])/math.sqrt( GR1[1]**2 + GR2[1]**2 )
        if ( drift_significance>n_sigma ):
            print 'Rejecting file: ', file_name
            print 'Significant time drifting detected'
            continue
        # save the result:
        benchmark_results[ results.identifier ] = results

# do other computations:
# sort the model keys:
sorted_benchmark_results_keys = sorted(benchmark_results)

# prune the benchmark results removing very similar results:
for ind, benchmark_name in enumerate(sorted_benchmark_results_keys):
    if ind==0:
        continue
    prev_key = sorted_benchmark_results_keys[ind-1]
    if benchmark_results[benchmark_name].is_close( benchmark_results[prev_key] ):
        del benchmark_results[prev_key]

# update sorted list:
sorted_benchmark_results_keys = sorted(benchmark_results)

# do derived calculations:
for res in sorted_benchmark_results_keys:
    benchmark_results[res].compute_relative_bench()
