clear;
clc;

addpath("C:\Users\김동욱\OneDrive\Desktop\GNSSLAB\2026-동계방학\Study_1\parsing\data");
addpath(genpath('data'))

[nav_data, nav_header] = rinexread("BRDC00IGS_R_20253590000_01D_MN.rnx");
[obs_data, obs_header] = rinexread("YONS00KOR_R_20253590000_01D_30S_MO.rnx");

%% Navigation Parsing %%
%======================%
%
%======================%


