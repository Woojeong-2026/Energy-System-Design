% MAIN_CHARGING Charging results, independent checks and regression tests.
% MATLAB R2026a + Cantera 3.0.0 toolbox. Run from the repository root.
% The full assignment entry point remains main.m (includes discharging/RTE).
clear; clc; close all;
here=fileparts(mfilename('fullpath'));
addpath(fullfile(here,'charging'),fullfile(here,'tests'));
chargingSummary=run_charging_review(here);
test_charging_review(here);
fprintf('Charging review complete: %s\n',fullfile(here,'results','charging_review'));
