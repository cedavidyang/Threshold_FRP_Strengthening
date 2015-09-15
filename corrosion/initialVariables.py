import numpy as np
from constants import *
# initiatin of variables for Likelihood Weighting Sampling

# initial parameters

# Structural age
service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL) 

weight_sum = 0
n_iter = 1
seed_indx = np.arange(1,33,2)
    
chloride_sums = np.zeros(1)
corrosion_state_sums = np.zeros(1)
corrosion_rate_sums = np.zeros(1)
mean_corrosion_rate_sums = np.zeros(1)
residual_diameter_sums = np.zeros(1)
radial_pressure_sums = np.zeros(1)
crack_prob_sums = np.zeros(1)
ds_crack_sums = np.zeros(1)
crack_width_sums = np.zeros(1)
diffusion_crack_sums = np.zeros(1)
    
chloride_data = np.array([]).reshape(service_time.size,0)
corrosion_state_data = np.array([]).reshape(service_time.size,0)
corrosion_rate_data = np.array([]).reshape(service_time.size,0)
mean_corrosion_rate_data = np.array([]).reshape(service_time.size,0)
residual_diameter_data = np.array([]).reshape(service_time.size,0)
radial_pressure_data = np.array([]).reshape(service_time.size,0)
crack_initiation_data = np.array([]).reshape(service_time.size,0)
crack_width_data = np.array([]).reshape(service_time.size,0)
diffusion_crack_data = np.array([]).reshape(service_time.size,0)
likelihood_weighting_data = np.array([])