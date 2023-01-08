%% MODEL COMPARISON
% script to plot model temperature output together with real data. 
% uses function mymodelcomp.m
% Information about netCDF-files: 
% Variable is temperature, T(lon, lat, time) [K] for different CO2 concentrations 
% longitudes go only in east direction 
% latitudes go from south (negative) to north (positive)
% grid resolution: 3x3 degrees
% time given in days
% Drilling site coordinates (based on Vetteroretti et al. (2022)
% NGRIP: 75.10N, 42.32W
% GRIP: 72.58N, 37.64W
% GISP2: 72.60N, 38.50W
% NEEM: 77.45, 51.06W
% WDC: 79.3S, 112W
% TAL: 72.5S, 159.0E 
% DOME F: 77.2S; 39.4E
% EDC: 75.0S; 123.1E
% EDML: 75.0S, 0.0E

clear variables; close all; clc
warning('off','all')
set(0,'DefaultFigureVisible','on');

% DEFINE COLOUR CODES
cNGRIP = [0 0.4470 0.7410];
cGRIP = [0.7350 0.0780 0.1840];
cNEEM = [0.9290 0.6940 0.1250];
cGISP = [0.4660 0.6740 0.1880];
cWDC = [0.4940 0.1840 0.5560];
cTAL = [.2 .6 .4];
cDF = [0.6235 0.8235 0.9882];
cEDC = [0.8902 0.4627 0.1961];
cEDML = [0.8902 0.6039 0.9608];

% Coordinates of drilling sites

coord = table('Size',[9, 3],'VariableTypes',["string" 'double' 'double'],'VariableNames',["Core" "Longitude" "Latitude"]);
coord(:,1) = {'NGRIP', 'GRIP', 'NEEM', 'GISP2', 'WDC', 'TAL', 'DOME F', 'EDC', 'EDML'}.';
coord(:,2) = array2table([360-42.32 360-37.64 360-51.06 360-38.50 360-112.086 159.06 39.40 123.1 0.0684]');
coord(:,3) = array2table([75.1 72.58 77.45 72.60 -79.3 -72.48 -77.2 -75.0 -75.0]');

%% 185 ppm 

file = 'cesmi6gat31rblc185_ANN_210912_998911_cam2_decclimots_TS.nc'; % load file
T = ncread(file,'TS'); 
lon = ncread(file,'lon');
lat = ncread(file,'lat'); 
time = ncread(file,'time')/365; % time in years

mymodelcomp(coord,T,lon,lat,time,'185 ppm', [3500 5000], [-51 -39]) %, stacks_comp);
%% 200 ppm

file = 'cesmi6gat31rblc200_ANN_210912_998911_cam2_decclimots_TS.nc'; % load file
T = ncread(file,'TS'); % Temperature = TS(lon, lat, time) [K]
lon = ncread(file,'lon'); % only EAST! from 0 to 360°
lat = ncread(file,'lat'); % from South (neg) to North (pos)
time_200 = ncread(file,'time')/365; % UNIT???? DAYS????
d18Otimed = readtable('d18OtimedGI-5.2.xlsx');

mymodelcomp(coord,T,lon,lat,time_200,'200 ppm',[4000 6000], [-50 -38]) %, stacks_comp, d18Otimed);

GLstack= zeros(length(d18Otimed.NGRIP),1);
for i = 1:length(d18Otimed.NGRIP)
    GLstack(i) = (d18Otimed.NGRIP(i) + d18Otimed.GRIP(i) + d18Otimed.NEEM(i))/3;
end 

%% 210 ppm

file = 'cesmi6gat31rblc210_ANN_210912_998911_cam2_decclimots_TS.nc'; 
T = ncread(file,'TS'); 
lon = ncread(file,'lon'); 
lat = ncread(file,'lat'); 
time = ncread(file,'time')/365; 

mymodelcomp(coord,T,lon,lat,time,'210 ppm',[4000 6000], [-51 -38]); %, stacks_comp);

%% HEINRICH 210 ppm

file = 'cesmi6gat31rblc210_Heinrich_ANN_210912_998911_cam2_decclimots_TS.nc'; 
T = ncread(file,'TS'); 
lon = ncread(file,'lon'); 
lat = ncread(file,'lat');
time = ncread(file,'time')/365; 
mymodelcomp(coord,T,lon,lat,time,'H210 ppm',[2500 5000], [-51 -38]); %, stacks_comp)
%% 220 ppm

file = 'cesmi6gat31rblc220_ANN_210912_998911_cam2_decclimots_TS.nc'; % load file
T = ncread(file,'TS'); % Temperature = TS(lon, lat, time) [K]
lon = ncread(file,'lon'); % only EAST! from 0 to 360°
lat = ncread(file,'lat'); % from South (neg) to North (pos)
time_220 = ncread(file,'time')/365; 
mymodelcomp(coord,T,lon,lat,time_220,'220 ppm',[5000 7000], [-51 -37]) %, stacks_comp);

%% 225 ppm

file = 'cesmi6gat31rblc225_ANN_210912_998911_cam2_decclimots_TS.nc'; % load file
T = ncread(file,'TS'); 
lon = ncread(file,'lon'); 
lat = ncread(file,'lat'); 
time = ncread(file,'time')/365; 
mymodelcomp(coord,T,lon,lat,time,'225 ppm',[5000 7500],[-51 -38]); %, stacks_comp);

%% 230 ppm

file = 'cesmi6gat31rblc230_ANN_210912_998911_cam2_decclimots_TS.nc'; 
T = ncread(file,'TS'); 
lon = ncread(file,'lon'); 
lat = ncread(file,'lat'); 
time = ncread(file,'time')/365; 
mymodelcomp(coord,T,lon,lat,time,'230 ppm',[time(1), 2600], [-51 -36]); %, stacks_comp);
