# Script to train and save stan models
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 14:46:24 2019

@author: Christof Naumzik
"""

import pystan
from pystan import StanModel 
import stan_utility
import os
import sys
import numpy as np
os.chdir('/local/home/cnaumzik/FakeNews')
model_name = sys.argv[1]
data_set = sys.argv[2]
states = sys.argv[3]
print('Loading' + data_set + 'now')
data = pystan.read_rdump('Data/'+ data_set +'_data_dump.R')
print('Dataset loaded')
print('Set states to ' + states)
data["S"] = np.int(states)
print("Fitting model " + model_name + " now")
model = StanModel('Code/StanCode/' + model_name + '.stan')

print("Fitting full model")
fit = model.sampling(data = data, init_r = 1, refresh = 10, iter = 2500,chains = 2)

print('Model fit successful - Saving fit now')

stan_utility.save_fit.stanfit_to_hdf5(fit = fit, file_name = "Data/Model fits/" + model_name +  "_" + states +  "_" + data_set)


