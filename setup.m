
% General setup for modeling
% Copyright 2021 Vibha Viswanathan. All rights reserved. 

% Model params
f_low = 125;
f_high = 8000;
nfilts = 30; % there are 30 abutting channels between 125 and 8000 (cams(8000)-cams(125))
list_CFs = invcams(linspace(cams(f_low), cams(f_high), nfilts)); % Equally space CFs on ERB scale
nCFs = numel(list_CFs);
fs_model = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz) 
cohc  = 1.0;    % normal ohc function
cihc  = 1.0;    % normal ihc function
% Species: 1 for cat (2 for human with Shera et al. tuning; 3 for human
% with Glasberg & Moore tuning)
species = 2;
noiseType = 0;  % 1 for variable fGn (0 for fixed fGn or const. spont rate)
spont = 10; % specify AN fiber spont rate 
% implnt: "0" for approximate or "1" for actual implementation of the
% power-law functions in the Synapse (i.e., neural adaptation)
implnt = 0;
tabs = 0.6e-3; % absolute refractory period
trel = 0.6e-3; % relative refractory period

% Stimulus params
nrep = 200; % number of stimulus repetitions (e.g., 50);
fs_STeVI = 48000; % sampling rate of speech stimuli
targetSPL = 60; % dB SPL of target speech, produces sufficient response
% for consonants in quiet, and also does not saturate the nerve response for
% the loudest stimulus (i.e., SiSNM or SiDCmod at -18 dB SNR).
consonants = {'b','ch','d','dh','f','g','j','k','l','m',...
    'n','p','r','s','sh','t','th','v','z','zh',};
nconson = numel(consonants);
conditions = {'SiQuiet','SiSSN_-8dB','SiB_-8dB','SiSNM_-18dB','vocodedSiB_0dB',...
    'vocodedSiQuiet'}; % Note: SiSNM is the same as SiDCmod
CondNamesNoSNR = {'SiQuiet','SiSSN','SiB','SiSnM','Vocoded SiB','Vocoded SiQuiet'};
nconds = numel(conditions);
talkers = {'Female1','Female2','Male1','Male2'};
ntalkers = numel(talkers);
load('consonantSegmentation.mat'); % Load segmentation info (variable startsec_all)
Cdur_sec = 0.1042; % consonant duration in seconds (5000 samples in fs_STeVI)
T = 2; % time duration of stimulus, model simulation, and resulting PSTH (seconds)
targRMS = 0.02; % Target speech starts out at RMS of 0.02.

% PSTH params
psthbinwidth = 1e-3; % binwidth in seconds;
nTimePtsPerPSTHbin = round(psthbinwidth*fs_model);  % number of time points per psth bin
numPSTHbins = round((T*fs_model)/nTimePtsPerPSTHbin); % number of PSTH bins total
fs_PSTH = 1/psthbinwidth;
consBinsStart = ceil(startsec_all./psthbinwidth); % in PSTH bins
Cdur_bins = ceil(Cdur_sec./psthbinwidth); % consonant duration (in PSTH bins)

% Paths
ANsimResultPath = '../results/modeling/';
perceptualDatapath = '../results/online_main_study';
modelPath = '../models/BEZ2018model/'; % path to Bruce et al. (2018) model  
addpath(modelPath);
stimPath = '../stimuli/'; % path to audio stimuli

% Define modulation filterbank per Relano-Iborra et al. (2016):
% third-order low-pass filter with a cutoff frequency fc=1Hz in parallel with eight
% second-order bandpass filters with octave spacing, a constant quality factor Q
% of 1, and center frequencies ranging from 2 to 256 Hz.
cf_mod = [1 2 4 8 16 32 64 128 256]; % center frequencies of the modulation filters
nModBands = numel(cf_mod);

