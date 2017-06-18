% Main function of the SPAC model
% Date, citation... 
clear;clc;close all;

% Model input, including time series of climate and LAI, 
% vegetation and soil properties, and initial conditions
load Input.mat; 

% Output path, stores hourly states and fluxes
outpath = 'Output\'; delete([outpath,'*.txt']);
% 
NumofDay = 365; % model time length in days
Ppt = -1; % Stochastic precipitation will be generated when Ppt == -1; 
          % Otherwise provide daily precipitation time series here. 
          
h = solveSPAC(outpath,CLM,Soil,Veg,IC,Const,NumofDay,dLAI,Ppt);

plotFigures(outpath,Veg) % Plot stomatal conductance and hydraulic architecture 