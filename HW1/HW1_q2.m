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
%% Analysis time window: [-0.2, 0.8] s

windowIdx = time >= -0.2 & time <= 0.8; % this creates a logic array (true/false mask)
timeWindow = time(windowIdx); % if the corresponding idx has true label, time value is kept.

fprintf('Analysis window: %.4f to %.4f s\n', ...
    timeWindow(1), timeWindow(end));

%% Grand average at Cz for each magnitude

magnitudes = [0 3 6 9 12]; % it's like the degree of rotation error applied to the user's cursor movement

grandAvg = zeros(length(timeWindow), length(magnitudes));

figure;
hold on;

for m = 1:length(magnitudes)
    trialIdx = magnitude == magnitudes(m); % create a true/false mask that trial's corresponding magnitude matches current one
    % time x trials
    % windowIdx -> select time from -0.2 to +0.8 s 
    % czIdx -> select only the Cz electrode

    czData = squeeze(rotation_data(windowIdx, czIdx, trialIdx));

    % Average across trials
    grandAvg(:,m) = mean(czData, 2); % dimension starts from 1, not 0

    plot(timeWindow, grandAvg(:,m), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d deg', magnitudes(m)));

    fprintf('%2d deg: %d trials\n', ...
        magnitudes(m), sum(trialIdx));
end

xline(0, '--k', 'Trigger', 'HandleVisibility', 'off');

xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('Grand Average ErrP at Cz');
legend('show');
grid on;
%% Grand average across ALL error trials at Cz

errorIdx = label == 1; % it choose error data (error trial has label 1)
errorCz = squeeze( ...
    rotation_data(windowIdx, czIdx, errorIdx));

% time x error trials -> time x 1
errorGrandCz = mean(errorCz, 2);

figure;
plot(timeWindow, errorGrandCz, ...
    'LineWidth', 2);

xline(0, '--k');

xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('Grand Average of All Error Trials at Cz');
grid on;

%% Identify ERN(first negative peak) and Pe(first positive peak)
postIdx = timeWindow >= 0;

% only use signal after trigger onset
postTime = timeWindow(postIdx);
postSignal = errorGrandCz(postIdx);

% findpeaks only finds positive peak, if we want to find negative peak,
% give it a neg sign~
% findpeaks find all local peak in the data
% I choose the first one!
[~, negLocs] = findpeaks(-postSignal); 
if isempty(negLocs)
    error('No negative peak found after trigger.');
end
ERN_idx_post = negLocs(1);
ERN_time = postTime(ERN_idx_post);
ERN_amp  = postSignal(ERN_idx_post);

[~, posLocs] = findpeaks(postSignal);

validPe = posLocs(posLocs > ERN_idx_post);

if isempty(validPe)
    error('No positive peak found after ERN.');
end

Pe_idx_post = validPe(1);

Pe_time = postTime(Pe_idx_post);
Pe_amp  = postSignal(Pe_idx_post);
fprintf('\nERN = %.1f ms, amplitude = %.3f uV\n', ...
    ERN_time*1000, ERN_amp);

fprintf('Pe  = %.1f ms, amplitude = %.3f uV\n', ...
    Pe_time*1000, Pe_amp);

figure;

% EEG waveform
plot(timeWindow, errorGrandCz, ...
    'LineWidth', 2.5, ...
    'Color', [0.2 0.7 1.0]);   % bright blue

hold on;

% ERN: bright red circle
plot(ERN_time, ERN_amp, 'o', ...
    'MarkerSize', 10, ...
    'LineWidth', 2, ...
    'MarkerEdgeColor', [1 0.3 0.3], ...
    'MarkerFaceColor', [1 0.3 0.3]);

% Pe: bright green square
plot(Pe_time, Pe_amp, 's', ...
    'MarkerSize', 10, ...
    'LineWidth', 2, ...
    'MarkerEdgeColor', [0.3 1.0 0.4], ...
    'MarkerFaceColor', [0.3 1.0 0.4]);

% Trigger
xline(0, '--', ...
    'Color', [1 1 1], ...
    'LineWidth', 1.2);

xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('ERN and Pe at Cz');

legend('Error Grand Average', ...
       'ERN', 'Pe', 'Trigger');

grid on;
%% Convert peak times to samples in original epoch
% min() returns [value, idx]
% ERN_sample is the index representing ERN happened
[~, ERN_sample] = min(abs(time - ERN_time));
[~, Pe_sample]  = min(abs(time - Pe_time));
fprintf('ERN sample = %d, time = %.4f s\n', ...
    ERN_sample, time(ERN_sample));

fprintf('Pe sample  = %d, time = %.4f s\n', ...
    Pe_sample, time(Pe_sample));

% channel x errorTrials
ERN_error = squeeze( ...
    rotation_data(ERN_sample, :, errorIdx));

% Average across error trials
ERN_topo = mean(ERN_error, 2);

Pe_error = squeeze( ...
    rotation_data(Pe_sample, :, errorIdx)); % get only error data and time exactly on Pe happened

Pe_topo = mean(Pe_error, 2);
% Both size(ERN_topo) & size(Pe_topo) is the size 32x1, cuz every electrode has the 1 exact same point as
% ERN/Pe happens

%% Topoplots before CAR

figure;

subplot(1,2,1);
topoplot(ERN_topo, chanlocs, 'maplimits', 'absmax');

title(sprintf('ERN - %.0f ms', ERN_time*1000), ...
    'Color', 'k', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

cb = colorbar;
cb.Color = 'k';
cb.FontSize = 10;

subplot(1,2,2);
topoplot(Pe_topo, chanlocs, 'maplimits', 'absmax');

title(sprintf('Pe - %.0f ms', Pe_time*1000), ...
    'Color', 'k', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

cb = colorbar;
cb.Color = 'k';
cb.FontSize = 10;
%% Perform Common Average Reference (CAR)
% data dimension is: 1024 × 32 × 320 , calculate avg of channels(dim 2) as
% CAR
rotation_data_car = ...
    rotation_data - mean(rotation_data, 2);
size(rotation_data_car)



%% Grand average at Cz after CAR

grandAvgCAR = zeros(length(timeWindow), length(magnitudes));

figure;
hold on;

for m = 1:length(magnitudes)

    trialIdx = magnitude == magnitudes(m);

    czDataCAR = squeeze( ...
        rotation_data_car(windowIdx, czIdx, trialIdx));

    grandAvgCAR(:,m) = mean(czDataCAR, 2);

    plot(timeWindow, grandAvgCAR(:,m), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d deg', magnitudes(m)));
end

xline(0, '--k', 'Trigger');

xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('Grand Average ErrP at Cz after CAR');
legend('show');
grid on;

%% CAR topoplots at same ERN / Pe time points

ERN_error_car = squeeze( ...
    rotation_data_car(ERN_sample, :, errorIdx));

Pe_error_car = squeeze( ...
    rotation_data_car(Pe_sample, :, errorIdx));

ERN_topo_car = mean(ERN_error_car, 2);
Pe_topo_car  = mean(Pe_error_car, 2);


figure;

%% ERN
subplot(1,2,1);

topoplot(ERN_topo_car, chanlocs, ...
    'maplimits', 'absmax');

title(sprintf('ERN after CAR - %.0f ms', ERN_time*1000), ...
    'Color', 'k', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

cb1 = colorbar;
cb1.Color = 'k';       % colorbar text/ticks = black
cb1.FontSize = 10;


%% Pe
subplot(1,2,2);

topoplot(Pe_topo_car, chanlocs, ...
    'maplimits', 'absmax');

title(sprintf('Pe after CAR - %.0f ms', Pe_time*1000), ...
    'Color', 'k', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

cb2 = colorbar;
cb2.Color = 'k';       % colorbar text/ticks = black
cb2.FontSize = 10;


% Figure background
set(gcf, 'Color', 'w');

%% CCA spatial filtering

% xdata has dimension: time × 32 channels × 320 trials
Xdata = rotation_data(windowIdx,:,:);

nTimeWindow = size(Xdata,1);
nChannels   = size(Xdata,2);
nTrials     = size(Xdata,3);

% Error template: time x channels
errorTemplate = mean(Xdata(:,:,label == 1), 3);

% Correct template: time x channels
correctTemplate = mean(Xdata(:,:,label == 0), 3);

X = zeros(nTimeWindow*nTrials, nChannels);
Y = zeros(nTimeWindow*nTrials, nChannels);

% stack every trial vertically 
for tr = 1:nTrials

    rows = (tr-1)*nTimeWindow + ...
        (1:nTimeWindow);

    % Actual EEG for this trial
    X(rows,:) = Xdata(:,:,tr);

    % Corresponding class template
    if label(tr) == 1
        Y(rows,:) = errorTemplate;
    else
        Y(rows,:) = correctTemplate;
    end
end
% size(X)
% size(Y)

[A, B, r, U, V] = canoncorr(X, Y);
% 
% size(A) = 32 x 32
fprintf('\nFirst five canonical correlations:\n');

for k = 1:5
    fprintf('CCA %d: %.4f\n', k, r(k));
end
%% CCA topoplots

figure;

for k = 1:5

    subplot(2,3,k);

    topoplot(A(:,k), chanlocs, ...
        'maplimits', 'absmax');

    title(sprintf('CCA Component %d', k), ...
        'Color','k', ...
        'FontWeight','bold');

    % Colorbar with black text
    cb = colorbar;
    cb.Color = 'k';
end

% Optional: white figure background
set(gcf, 'Color', 'w');
%% topoplot
% topoplot(data vector,selectedChannels) % data vector [n of channels x 1]


%% grand average plot
% 
% %% CCA
% % maximize correlation between the single-trial data and the class grand average template
% 
% [spatialFilter,~] =  canoncorr(concat_data', concat_ga'); % spatial filter [n channels x n components]

