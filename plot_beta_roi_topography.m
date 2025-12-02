function plot_beta_roi_topography()
% PLOT_BETA_ROI_TOPOGRAPHY
%
% This script computes and plots scalp topographies of beta-band ERD, ERS,
% and peak-to-peak (ERS - ERD) for broad frontal, left and right ROIs
% for selected subjects and conditions.
%
% It was used as an exploratory tool to identify channels showing maximal
% movement-related beta modulation and to guide the definition of ROIs in:
%
%   Tatti et al., "Insights into Central Fatigue in Multiple Sclerosis:
%   Movement-Related Beta Oscillatory Activity and its Association with
%   Brain Neurophysiology and Structure" (Brain Communications).
%
% REQUIREMENTS
%   - MATLAB R20xx
%   - FieldTrip installed and added to the MATLAB path
%   - EEGLAB (or at least EEGLAB's topoplot.m) on the MATLAB path
%   - Preprocessed EEG data in FieldTrip format:
%       data.elec  : electrode positions
%       data.trial : trials x channels x time
%       data.time  : 1 x time vector
%   - chanlocs.mat file containing EEGLAB channel locations (variable chanlocs)
%
% The script:
%   1) Computes a wavelet time–frequency analysis in the 13–25 Hz beta range.
%   2) Averages beta power across frequency -> trials x channels x time.
%   3) Uses three broad ROIs (frontal, left, right) to find ERD and ERS
%      peak times per trial.
%   4) Normalizes ERD/ERS power by a mean beta reference for each channel.
%   5) Plots topographies for ERD, ERS, and peak-to-peak for each ROI.

%% ====================== USER SETTINGS ==================================

% Root folder of your project (adapt!)
projectRoot = '/path/to/TAPPING_EEG_MS_Siena/TAPPING_EEG';

% Folder containing preprocessed FieldTrip data
% expected: <projectRoot>/<cond>/Preprocessing_avgref_ALE/Fieldtrip/Sbj<s>_<cond>_fixchan.mat
dataRoot = projectRoot;

% chanlocs file (EEGLAB format) with variable 'chanlocs'
chanlocsFile = fullfile(projectRoot, 'Scripts', 'chanlocs.mat');  % <-- change if needed

% Broad ROIs (indices refer to your chanlocs / FieldTrip channel order)
% These are the "broad" zones you used to find typical ERD/ERS latencies.
ROI_frontal = [10 11 12 13 9 8 7 5 4 61 19 20 21 22 18 17];
ROI_left    = [28 29 30 31 40 39 38 37 49 48 47 46];
ROI_right   = [28 27 26 25 37 36 35 34 46 45 44 43];

% Condition and subjects to explore
% (change these for your study; here is just one example)
cond     = '3';        % '1' = MS_F, '2' = MS_NF, '3' = HC (for example)
groupTag = 'HC';       % only used for naming figures, not strictly necessary
subjects = 5;          % e.g. [1 2 3 4] or just 5 as in your original example

% Time–frequency analysis settings
f_o_i = 13:0.5:25;     % beta band frequencies (Hz)
t_o_i = -0.5:0.01:1.5; % time of interest (s)
cycles = 3:10;         % used to compute width in the original code

% ERD and ERS search windows in indices (relative to t_o_i discretization)
% These are taken from your original script:
erdTimeIdx   = 51:121;    % "early" beta suppression window
ersTimeIdx   = 121:151;   % "late" beta rebound window

%% ====================== PREPARATION ====================================

% Load chanlocs for topoplot
if ~exist(chanlocsFile, 'file')
    error('chanlocs file not found: %s', chanlocsFile);
end
load(chanlocsFile, 'chanlocs');

% Make sure EEGLAB / topoplot is on the path
% (you can uncomment eeglab if you want it to launch automatically)
% eeglab; close all;

count = 1; % subject counter if you later want to store ROIs interactively

%% ====================== LOOP OVER SUBJECTS =============================

for s = subjects

    fprintf('\n=== Subject %d, condition %s (%s) ===\n', s, cond, groupTag);

    % Load preprocessed FieldTrip data
    dataFile = fullfile(dataRoot, cond, 'Preprocessing_avgref_ALE', 'Fieldtrip', ...
                        sprintf('Sbj%d_%s_fixchan.mat', s, cond));
    if ~exist(dataFile, 'file')
        warning('File not found: %s. Skipping subject.', dataFile);
        continue;
    end
    load(dataFile, 'data'); % loads variable 'data'

    %% Prepare layout just in case (not strictly required for this script)
    cfg            = [];
    cfg.elec       = data.elec;
    cfg.rotate     = [];
    cfg.projection = 'orthographic';
    layout         = ft_prepare_layout(cfg, data); %#ok<NASGU>

    %% Time–frequency analysis in the 13–25 Hz beta range

    cfg            = [];
    cfg.method     = 'wavelet';
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.trials     = 'all';
    cfg.keeptrials = 'yes';
    cfg.keeptapers = 'no';
    cfg.foi        = f_o_i;
    cfg.toi        = t_o_i;
    cfg.width      = ceil(f_o_i / (numel(cycles) - 1)) + 2; % original formula

    freq = ft_freqanalysis(cfg, data);

    % Average over beta frequencies: trials x channels x time
    % freq.powspctrm: trials x channels x freqs x time
    meanbetafreq = squeeze(mean(freq.powspctrm, 3));

    %% Compute ROI-mean time courses and ERD/ERS peak times

    roi1mean = squeeze(nanmean(meanbetafreq(:, ROI_frontal, :), 2)); % frontal
    roi2mean = squeeze(nanmean(meanbetafreq(:, ROI_left,    :), 2)); % left
    roi3mean = squeeze(nanmean(meanbetafreq(:, ROI_right,   :), 2)); % right

    % For each trial, find ERD (min) and ERS (max) within the specified windows
    [~, edT1] = min(roi1mean(:, erdTimeIdx), [], 2);
    [~, esT1] = max(roi1mean(:, ersTimeIdx), [], 2);

    [~, edT2] = min(roi2mean(:, erdTimeIdx), [], 2);
    [~, esT2] = max(roi2mean(:, ersTimeIdx), [], 2);

    [~, edT3] = min(roi3mean(:, erdTimeIdx), [], 2);
    [~, esT3] = max(roi3mean(:, ersTimeIdx), [], 2);

    % Convert local indices to full time indices (as in the original script)
    erdT1 = edT1 + erdTimeIdx(1) - 1;
    ersT1 = esT1 + ersTimeIdx(1) - 1;
    erdT2 = edT2 + erdTimeIdx(1) - 1;
    ersT2 = esT2 + ersTimeIdx(1) - 1;
    erdT3 = edT3 + erdTimeIdx(1) - 1;
    ersT3 = esT3 + ersTimeIdx(1) - 1;

    % Reference: mean over all time and trials for each channel
    % ref: 1 x nChannels
    ref = nanmean(nanmean(meanbetafreq, 1), 3);

    nTrials   = size(meanbetafreq, 1);
    nChannels = size(meanbetafreq, 2);

    % Preallocate matrices for normalized ERD/ERS per trial and channel
    erD1sTotal = nan(nTrials, nChannels);
    erS1sTotal = nan(nTrials, nChannels);
    erD2sTotal = nan(nTrials, nChannels);
    erS2sTotal = nan(nTrials, nChannels);
    erD3sTotal = nan(nTrials, nChannels);
    erS3sTotal = nan(nTrials, nChannels);

    % Normalization over time: value at ERD/ERS peak divided by channel reference
    for ch = 1:nChannels
        for tr = 1:nTrials
            erD1sTotal(tr, ch) = meanbetafreq(tr, ch, erdT1(tr)) / ref(ch);
            erS1sTotal(tr, ch) = meanbetafreq(tr, ch, ersT1(tr)) / ref(ch);

            erD2sTotal(tr, ch) = meanbetafreq(tr, ch, erdT2(tr)) / ref(ch);
            erS2sTotal(tr, ch) = meanbetafreq(tr, ch, ersT2(tr)) / ref(ch);

            erD3sTotal(tr, ch) = meanbetafreq(tr, ch, erdT3(tr)) / ref(ch);
            erS3sTotal(tr, ch) = meanbetafreq(tr, ch, ersT3(tr)) / ref(ch);
        end
    end

    % Average across trials
    erS1s = mean(erS1sTotal, 1);
    erD1s = mean(erD1sTotal, 1);
    p2p1  = mean(erS1sTotal - erD1sTotal, 1);

    erS2s = mean(erS2sTotal, 1);
    erD2s = mean(erD2sTotal, 1);
    p2p2  = mean(erS2sTotal - erD2sTotal, 1);

    erS3s = mean(erS3sTotal, 1);
    erD3s = mean(erD3sTotal,  1);
    p2p3  = mean(erS3sTotal - erD3sTotal, 1);

    % Map limits for plots
    maxS1 = max(erS1s);
    minD1 = min(erD1s);
    maxP1 = max(p2p1);

    maxS2 = max(erS2s);
    minD2 = min(erD2s);
    maxP2 = max(p2p2);

    maxS3 = max(erS3s);
    minD3 = min(erD3s);
    maxP3 = max(p2p3);

    %% Plot topographies for each ROI (ERD, ERS, P2P)
    % For figure naming it is convenient to encode subject and group

    figure('Name', sprintf('Sub%d_%s_ROI1_ERD', s, groupTag));
    topoplot(erD1s, chanlocs, 'maplimits', [minD1 1], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Frontal ROI ERD', s, groupTag));

    figure('Name', sprintf('Sub%d_%s_ROI1_ERS', s, groupTag));
    topoplot(erS1s, chanlocs, 'maplimits', [1 maxS1], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Frontal ROI ERS', s, groupTag));

    figure('Name', sprintf('Sub%d_%s_ROI1_P2P', s, groupTag));
    topoplot(p2p1, chanlocs, 'maplimits', [0 maxP1], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Frontal ROI P2P', s, groupTag));

    % Left ROI
    figure('Name', sprintf('Sub%d_%s_ROI2_ERD', s, groupTag));
    topoplot(erD2s, chanlocs, 'maplimits', [minD2 1], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Left ROI ERD', s, groupTag));

    figure('Name', sprintf('Sub%d_%s_ROI2_ERS', s, groupTag));
    topoplot(erS2s, chanlocs, 'maplimits', [1 maxS2], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Left ROI ERS', s, groupTag));

    figure('Name', sprintf('Sub%d_%s_ROI2_P2P', s, groupTag));
    topoplot(p2p2, chanlocs, 'maplimits', [0 maxP2], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Left ROI P2P', s, groupTag));

    % Right ROI
    figure('Name', sprintf('Sub%d_%s_ROI3_ERD', s, groupTag));
    topoplot(erD3s, chanlocs, 'maplimits', [minD3 1], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Right ROI ERD', s, groupTag));

    figure('Name', sprintf('Sub%d_%s_ROI3_ERS', s, groupTag));
    topoplot(erS3s, chanlocs, 'maplimits', [1 maxS3], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Right ROI ERS', s, groupTag));

    figure('Name', sprintf('Sub%d_%s_ROI3_P2P', s, groupTag));
    topoplot(p2p3, chanlocs, 'maplimits', [0 maxP3], ...
             'gridscale', 100, 'headrad', 0.61, 'electrodes', 'numbers');
    title(sprintf('Subject %d %s - Right ROI P2P', s, groupTag));

    fprintf('Finished subject %d.\n', s);

    % OPTIONAL: interactive ROI selection (if you want to recreate how you picked channels)
    %{
    % After inspecting the topoplots, you can manually type the list of
    % channel indices for each ROI and store them in ROI_dx, ROI_sx, ROI_front.

    fprintf('Subject %d: see topoplots and type channel indices if you want to store new ROIs.\n', s);

    xy  = input('Enter channel indices for RIGHT ROI (e.g. [28 27 26 ...]): ', 's');
    ROI_dx(count, :) = str2num(xy); %#ok<ST2NM>

    xy  = input('Enter channel indices for LEFT ROI: ', 's');
    ROI_sx(count, :) = str2num(xy);

    xy  = input('Enter channel indices for FRONTAL ROI: ', 's');
    ROI_front(count, :) = str2num(xy);

    count = count + 1;
    %}

end % subjects

end
