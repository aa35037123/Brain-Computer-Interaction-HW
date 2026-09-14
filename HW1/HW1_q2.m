% [eeg,header] = sload('Subject_006_TESS_Online__feedback__s006_r001_2021_08_30_163523.gdf');
clear;
clc;
close all;
%% Install/setup BioSig

original_dir = pwd;   % remember current directory

cd('..');             % move up one directory
biosig_installer;       % or biosig_installer, depending on filename

cd(original_dir);     % come back

% Add folder containing sload.m
addpath(['ErrP_data_scripts'])
addpath(['MI_data_scripts']);
% Setup and load data
load('ErrP_data_scripts/ErrP_data_HW1.mat');
load('ErrP_data_scripts/ErrP_channels.mat');
% Extract data
rotation_data = trainingEpochs.rotation_data;   % time x channel x trial
label         = trainingEpochs.label;           % 320 x 1
magnitude     = trainingEpochs.magnitude;       % 320 x 1
fileID        = trainingEpochs.fileID;
sessionID     = trainingEpochs.sessionID;

% Dataset parameters
fs       = params.fsamp;        % 512 Hz
time     = params.epochTime;    % time vector
chanlocs = params.chanlocs;     % 1 x 32 channel-location struct
czIdx    = params.channelPlot;  % 15

fprintf('Data size: %d time x %d channels x %d trials\n', ...
    size(rotation_data,1), ...
    size(rotation_data,2), ...
    size(rotation_data,3));

fprintf('Sampling rate: %d Hz\n', fs);
fprintf('Cz channel index: %d\n', czIdx);

disp('Labels:');
disp(unique(label)');

disp('Magnitudes:');
disp(unique(magnitude)');
%% topoplot
% topoplot(data vector,selectedChannels) % data vector [n of channels x 1]


%% grand average plot
% 
% %% CCA
% % maximize correlation between the single-trial data and the class grand average template
% 
% [spatialFilter,~] =  canoncorr(concat_data', concat_ga'); % spatial filter [n channels x n components]

