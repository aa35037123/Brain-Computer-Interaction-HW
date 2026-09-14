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
addpath(genpath('MI_data_scripts'));
files = dir('MI_data_scripts/Subject_006_Session_006_TESS_Online_Visual/Subject_006_TESS_Online__feedback__s006_r*.gdf');

all_RH_raw = {};
all_LH_raw = {};
all_RH_car = {};
all_LH_car = {};
RH_trials_total = 0;
LH_trials_total = 0;

for r = 1:length(files)
    filename = fullfile(files(r).folder, files(r).name);
    fprintf('Processing %s\n', files(r).name);
    [eeg, header] = sload(filename);
    % Events
    eventType = header.EVENT.TYP;
    eventPos = header.EVENT.POS;
    trialStart = eventPos(eventType == 1000);
    fprintf('Trail start at samples');
    fprintf('%d ',trialStart);
    fprintf('sample\n');
    % first 32 column are EEG
    eeg = eeg(:,1:32);
    % sampling frequency is 512 hz
    fs = header.SampleRate;
    
    %% temporal filtering
    order = 4;
    low = 8;
    high = 13;
    nyquist = fs/2;
    norm_cutoff_freq = [low high]/nyquist;
    [b,a] = butter(order,norm_cutoff_freq,'bandpass');
    eeg_mu = filtfilt(b, a, eeg);
    %% Spatial filtering - CAR
    common_avg = mean(eeg_mu, 2); % 2 means take average on second dimension(column)
    eeg_mu_car = eeg_mu - common_avg;
    
    %% Deduce the number of trial
    n_RH = sum(header.EVENT.TYP == 769);
    n_LH = sum(header.EVENT.TYP == 770);
    RH_trials_total = RH_trials_total + n_RH;
    LH_trials_total = LH_trials_total + n_LH;
    fprintf('Run %d: \n', r);
    fprintf('RH trials: %d, LH trails: %d\n', n_RH, n_LH);
    % fprintf('Total trials: %d\n', n_RH + n_LH);

    %% Trial Extraction
    % each of RH_raw/LH_raw is a last 0.5s of a trial 
    [RH_raw,RH_car] = extract_last_half_sec( ...
        eventType,eventPos,eeg_mu,eeg_mu_car,...
        7691,[7692 7693],fs);
    [LH_raw,LH_car] = extract_last_half_sec( ...
        eventType,eventPos,eeg_mu,eeg_mu_car,...
        7701,[7702 7703],fs);
    %% Combine runs
    all_RH_raw = [all_RH_raw RH_raw];
    all_LH_raw = [all_LH_raw LH_raw];

    all_RH_car = [all_RH_car RH_car];
    all_LH_car = [all_LH_car LH_car];

end

%% Number of trials
fprintf('Total number of Right Hand trials: %d\n', RH_trials_total);
fprintf('Total number of Left Hand trials: %d\n', LH_trials_total);
fprintf('Total number of trials: %d\n', RH_trials_total + LH_trials_total);
% %% spatial filtering
% 
% CAR or Laplacian
% 
% %% trial extraction
% 
% data [time samples x channels x trials]
% 
%% mu power 

% instantenous power = time samples ^ 2
% for each trial
%     average power = sum of instanteneous power/n of samples

RH_power_raw = compute_mu_power(all_RH_raw);
LH_power_raw = compute_mu_power(all_LH_raw);
RH_grand_raw = mean(RH_power_raw,1);
LH_grand_raw = mean(LH_power_raw,1);

RH_power_car = compute_mu_power(all_RH_car);
LH_power_car = compute_mu_power(all_LH_car);

RH_grand_car = mean(RH_power_car,1);
LH_grand_car = mean(LH_power_car,1);
% min(RH_grand_raw)
% max(RH_grand_raw)
% 
% min(RH_grand_car)
% max(RH_grand_car)
% 
% min(LH_grand_raw)
% max(LH_grand_raw)
% 
% min(LH_grand_car)
% max(LH_grand_car)

%% topoplot
% topoplot(data vector,selectedChannels) % data vector [n of channels x 1]
load('MI_data_scripts/selectedChannels.mat');
figure;

subplot(2,2,1);

topoplot(RH_grand_raw(:), selectedChannels);
title('RH - No Spatial Filter', ...
    'Color','k','FontWeight','bold','FontSize',12);
colorbar;

cb = colorbar;       
cb.Color = 'k';      % set color bar's color as black
cb.FontSize = 10;

subplot(2,2,2);
topoplot(RH_grand_car(:), selectedChannels);
title('RH - CAR', ...
    'Color','k','FontWeight','bold','FontSize',12);
colorbar;
cb = colorbar;       
cb.Color = 'k';      % set color bar's color as black
cb.FontSize = 10;

subplot(2,2,3);
topoplot(LH_grand_raw(:), selectedChannels);
title('LH - No Spatial Filter', ...
    'Color','k','FontWeight','bold','FontSize',12);
colorbar;
cb = colorbar;       
cb.Color = 'k';      % set color bar's color as black
cb.FontSize = 10;

subplot(2,2,4);
topoplot(LH_grand_car(:), selectedChannels);
title('LH - CAR', ...
    'Color','k','FontWeight','bold','FontSize',12);


colorbar;
cb = colorbar;       
cb.Color = 'k';      % set color bar's color as black
cb.FontSize = 10;

%% grand average plot
% 
% %% CCA
% % maximize correlation between the single-trial data and the class grand average template
% 
% [spatialFilter,~] =  canoncorr(concat_data', concat_ga'); % spatial filter [n channels x n components]


%% === Functions ===
function [epochs_raw, epochs_car] = extract_last_half_sec( ...
    eventType,eventPos,rawData,carData,...
    startTrigger,endTriggers,fs)
    startIdx = find(eventType == startTrigger); % find return idx instead of value
    
    epochs_raw = {};
    epochs_car = {};
    
    nSamples = round(0.5*fs);
    for k = 1:length(startIdx)
        idx = startIdx(k);
        nextEnd = find( ...
            ismember(eventType(idx+1:end),endTriggers),...
            1,'first');
        if isempty(nextEnd)
            continue
        end
        endEventIdx = idx + nextEnd;
        taskEnd = eventPos(endEventIdx);
        epochStart = taskEnd-nSamples+1;
        epochEnd = taskEnd;
        if epochStart >= 1 && epochEnd <= size(rawData,1)

            epochs_raw{end+1} = ...
                rawData(epochStart:epochEnd,:);

            epochs_car{end+1} = ...
                carData(epochStart:epochEnd,:);
        end

    end
end

% Compute the average of every sample in a trial, and gather trial data together 
function power = compute_mu_power(trials)

    nTrials = length(trials);
    nChannels = size(trials{1},2);
    
    power = zeros(nTrials,nChannels);

    for k = 1:nTrials
        power(k,:) = mean(trials{k}.^2,1); % Takes average of all samples inside 1 trial (sample, channel)
    end
end