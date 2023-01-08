This collection of MATLAB code is part of the Master's thesis "Bipolar Phasing of 
Abrupt Climate Change During the Last Ice Age" by Anna Maria Klüssendorf


The code performs the dertermination of the onsets of 25 Dansgaard-Oeschger (DO) Events
in the NGRIP ice core as described in that work. Based on (bipolar) volcanic synchronisation, 
the respective timing of the DO events is determined in three other Greenlandic and five Antarctic ice cores. 
It further performs the data stacking and calculates the bipolat 
phasing between Northern Hemisphere abrupt warming and Southern Hemispehere cooling. 
Finally, it determines the bipolar phasing in temperature simulations of 
the CCSM4 model.

%% DETERMINE ONSETS OF DO EVENTS
% AMK_DOonsets_NGRIP.m
% Current script to find transitions, compare results to Buizert et al. and
% INTIMATE, Capron et al. 2020 and Myrvoll-Nielsen et al. 2022 as well as stacking of transitions. Uses a NEW DATASET including the entire core.
% uses myfindtransitions_single.m and NGRIP isotope data.

%% BIPOLAR PHASING 
% AMK_bipolar_phasing.m
% Analysis of the climate signal across the transitions associated with Dansgaard-Oeschger events by stacking the events aligned at the onset and determining the Antarctic response time for all events, those occurring in the early and the late ice age, respectively, and minor as well as major events. 
% This script uses the following functions:
% - mystacking.m to stack the events for the individual ice core records
% - myclasscomparison.m to compare stacks of the different classifications
% This script uses d18O data from NGRIP, GRIP, GISP2, NEEM, WDC, Talos Dome, EDC, Dome F, EDML
% + 
% TRANSITIONS determined in AMK_DOonsets-NGRIP.m
% BIPOLAR MATCHPOINTS as in Svensson et al. 2020
% GICC05 timescale d18O data as in Rasmussen et al. 2014
% BIPOLAR VOLCANOES as in Svensson et al. 2020

%% MODEL
% AMK_modelcomp.m
% script to plot model temperature output together with real data
% and determine bipolar phasing in simulated temperature at approximate location
% of all considered drilling sites under different CO2 forcing. 
% uses function mymodelcomp.m
