%% COMPUTE_TFR_TAPPING_MS
%
% Time–frequency analysis for the MS tapping study (FieldTrip).
%
% Steps:
%   1) Loop over three groups:
%        cond=1 : MS_F   (fatigued)
%        cond=2 : MS_NF  (non-fatigued)
%        cond=3 : HC     (healthy controls)
%   2) For each subject:
%        - Load preprocessed FieldTrip data with corrected channel order
%          (SbjX_cond_fixchan.mat)
%        - Apply surface Laplacian
%        - Compute single-trial TFR with Morlet wavelets 1–90 Hz, -0.5–1.5 s
%        - Normalize power relative to the average over the whole epoch
%        - Average across trials (chan × freq × time)
%   3) Store each subject TFR as SbjX_Y (X=subject, Y=cond)
%   4) Compute group grand-averages (ft_freqgrandaverage)
%   5) Plot example single-subject and group TFRs (imagesc)
%
% CHANNEL ORDER CORRECTION (IMPORTANT):
%
%   Due to a cabling error between the EEG cap and the amplifier,
%   the recorded channel order was incorrect in a subset of participants
%   (oldchanlocs, defined in the fix script).
%
%   For those subjects, we used a dedicated preprocessing script
%   (e.g., "fixchanlocs_tapping_ms.m", shown in the repository) that:
%     - remapped each row of data.trial{t} to the intended 61-channel
%       10–20 montage;
%     - reordered data.elec.elecpos, data.elec.pnt, and data.elec.label
%       accordingly;
%     - rewrote data.label to match the canonical electrode order.
%
%   That script was applied only to subjects listed in "oldchanlocs" for
%   each condition, and then saved as:
%       Sbj<subject>_<cond>_fixchan.mat
%
%   For subjects without the cabling issue, the same script simply copied
%   the original data into data_new (with canonical labels) and saved it,
%   so that ALL subjects now have a *_fixchan.mat file.
%
%   This TFR script therefore assumes that channel order and positions are
%   already correct in SbjX_cond_fixchan.mat for ALL subjects and does not
%   perform any further reordering.
%
% This script was used in:
%   Tatti et al., "Insights into Central Fatigue in Multiple Sclerosis:
%   Movement-Related Beta Oscillatory Activity and its Association with
%   Brain Neurophysiology and Structure" (Brain Communications).
%
% REQUIREMENTS:
%   - MATLAB
%   - FieldTrip
%   - Layout.mat: FieldTrip layout structure (variable "layout")
%   - labels_uppercase.mat: cell array "label" with channel labels
%
% USAGE:
%   - Check/adjust paths and subject lists below.
%   - Run from MATLAB:  compute_tfr_tapping_ms

clear; clc;

%% User settings ---------------------------------------------------------

% Base folder containing the FieldTrip-preprocessed data
baseFolder = '/Users/leonardo/Dropbox (City College)/Tapping_EEG_MS_Siena/TAPPING_EEG/';

% Layout and labels (FieldTrip-style)
layoutFile = '/Users/leonardo/Dropbox (City College)/Tapping_EEG_MS_Siena/chanlocs/Layout.mat';
labelsFile = '/Users/leonardo/Dropbox (City College)/Tapping_EEG_MS_Siena/chanlocs/labels_uppercase.mat';

% Group definitions (cond: 1 = MS_F, 2 = MS_NF, 3 = HC)
group(1).cond     = '1';
group(1).name     = 'MS_F';
group(1).subjects = [1:9 11:19];          % 18 subjects; Sbj10 no EEG

group(2).cond     = '2';
group(2).name     = 'MS_NF';
group(2).subjects = [1 3:14 16:22];       % 20 subjects; Sbj15 no EEG; Sbj2 excluded

group(3).cond     = '3';
group(3).name     = 'HC';
group(3).subjects = [1:7 9:18];           % 17 subjects; Sbj8 excluded (wrong sampling rate)

% TFR settings
f_o_i = 1:0.05:90;       % frequencies of interest (Hz)
t_o_i = -0.5:0.01:1.5;   % time of interest (s)
freq_resolution_factor = 4;

%% Load layout and labels -----------------------------------------------

if ~exist(layoutFile, 'file')
    error('Layout file not found: %s', layoutFile);
end
load(layoutFile, 'layout');   % must define "layout"

if ~exist(labelsFile, 'file')
    error('Labels file not found: %s', labelsFile);
end
load(labelsFile, 'label');    % cell array of channel labels

%% Loop over groups and subjects ----------------------------------------

for c = 1:numel(group)

    cond       = group(c).cond;
    subjList   = group(c).subjects;
    condName   = group(c).name;

    folder_field = fullfile(baseFolder, cond, 'Preprocessing_avgref_ALE', 'Fieldtrip');
    outputTFR    = fullfile(folder_field, 'TFR');
    outputFig    = fullfile(folder_field, 'TFR_plot');
    if ~exist(outputTFR, 'dir'), mkdir(outputTFR); end
    if ~exist(outputFig, 'dir'), mkdir(outputFig); end

    fprintf('\n=== Condition %s (%s) ===\n', cond, condName);

    for s = subjList

        % IMPORTANT: this file is the output of the fixchanlocs script
        dataFile = fullfile(folder_field, sprintf('Sbj%d_%s_fixchan.mat', s, cond));
        if ~exist(dataFile, 'file')
            warning('File not found: %s (skipping subject %d, cond %s)', dataFile, s, cond);
            continue;
        end

        fprintf('Subject %d, cond %s: loading %s\n', s, cond, dataFile);
        load(dataFile, 'data');   % FieldTrip data structure with corrected channels

        % 1) Laplacian filter (surface Laplacian)
        cfg           = [];
        cfg.method    = 'laplacian';
        cfg.layout    = layout;
        cfg.feedback  = 'no';
        data_lap      = ft_preprocessing(cfg, data);

        % 2) Time–frequency analysis (single trials)
        cfg           = [];
        cfg.method    = 'wavelet';
        cfg.output    = 'pow';
        cfg.channel   = 'all';
        cfg.trials    = 'all';
        cfg.keeptrials= 'yes';
        cfg.keeptapers= 'no';
        cfg.foi       = f_o_i;
        cfg.toi       = t_o_i;

        cycles        = 3:10;
        cfg.width     = ceil(f_o_i / (numel(cycles) - 1)) * freq_resolution_factor + 2;

        % Use Laplacian-filtered data (as in the final analysis)
        freq          = ft_freqanalysis(cfg, data_lap);

        % 3) Baseline normalization across the whole epoch
        %    refblock: chan × freq, averaged over trials and time
        refblock     = squeeze(nanmean(nanmean(freq.powspctrm,1),4));

        % Replicate baseline across trials and time
        nTime        = size(freq.powspctrm,4);
        refblock_toi = repmat(refblock, [1 1 nTime]);

        dims         = size(freq.powspctrm);
        normtrials   = nan(dims(1), dims(2), dims(3), dims(4));
        nTrials      = dims(1);
        for t = 1:nTrials
            normtrials(t,:,:,:) = (squeeze(freq.powspctrm(t,:,:,:)) - refblock_toi) ./ refblock_toi;
        end

        freq_norm            = freq;
        freq_norm.powspctrm  = normtrials;

        % 4) Quick check plot for one subject (optional)
        cfg_plot             = [];
        cfg_plot.parameter   = 'powspctrm';
        cfg_plot.channel     = 'all';
        cfg_plot.xlim        = [-0.5 1.5];
        cfg_plot.ylim        = [1 90];
        figure('Name', sprintf('TFR Sbj%d cond%s', s, cond));
        ft_singleplotTFR(cfg_plot, freq_norm);
        saveas(gcf, fullfile(outputFig, sprintf('TFR_Sbj%d_cond%s.png', s, cond)));

        % 5) Average across trials -> chan × freq × time
        freq_norm_av           = freq_norm;
        freq_norm_av.powspctrm = squeeze(nanmean(freq_norm.powspctrm,1));
        freq_norm_av.dimord    = 'chan_freq_time';
        freq_norm_av.label     = label;

        % 6) Store subject structure in a variable named SbjX_Y, as in the original script
        eval(sprintf('Sbj%d_%s = freq_norm_av;', s, cond));

        % 7) Save TFR structure
        save(fullfile(outputTFR, sprintf('TFR_Sbj%d_cond%s.mat', s, cond)), ...
             sprintf('Sbj%d_%s', s, cond), '-v7.3');

        clear freq freq_norm freq_norm_av data data_lap normtrials refblock refblock_toi;

    end
end

%% Grand averages per group ----------------------------------------------

cfgGA                = [];
cfgGA.keepindividual = 'yes';
cfgGA.foilim         = [1 90];
cfgGA.channel        = 'all';
cfgGA.parameter      = 'powspctrm';

[grandavg_group1] = ft_freqgrandaverage(cfgGA , ...
    Sbj1_1 , Sbj2_1 , Sbj3_1 , Sbj4_1 , Sbj5_1 , Sbj6_1 , ...
    Sbj7_1 , Sbj8_1 , Sbj9_1 , Sbj11_1 , Sbj12_1 , Sbj13_1 , ...
    Sbj14_1 , Sbj15_1 , Sbj16_1 , Sbj17_1 , Sbj18_1 , Sbj19_1);

[grandavg_group2] = ft_freqgrandaverage(cfgGA , ...
    Sbj1_2 , Sbj3_2 , Sbj4_2 , Sbj5_2 , Sbj6_2 , Sbj7_2 , ...
    Sbj8_2 , Sbj9_2 , Sbj10_2, Sbj11_2 , Sbj12_2 , Sbj13_2 , ...
    Sbj14_2 , Sbj16_2 , Sbj17_2 , Sbj18_2 , Sbj19_2 , Sbj20_2 , ...
    Sbj21_2 , Sbj22_2);

[grandavg_group3] = ft_freqgrandaverage(cfgGA , ...
    Sbj1_3 , Sbj2_3 , Sbj3_3 , Sbj4_3 , Sbj5_3 , Sbj6_3 , ...
    Sbj7_3 , Sbj9_3 , Sbj10_3, Sbj11_3 , Sbj12_3 , Sbj13_3 , ...
    Sbj14_3 , Sbj15_3 , Sbj16_3 , Sbj17_3 , Sbj18_3);

%% Simple TFR maps for one subject and three group grand-averages --------

% Example subject (here Sbj5_1)
plot_sbj = squeeze(nanmean(Sbj5_1.powspctrm(:,7:179,36:171),1)); % freq × time
figure('Name','Single-subject TFR (Sbj5_1)');
imagesc(Sbj5_1.time(36:171), Sbj5_1.freq(7:179), plot_sbj);
set(gca,'YDir','normal');
colormap(jet);
ylim([4 90]);
xlim([-0.15 1.2]);
caxis([-0.5 0.5]);
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Sbj5\_1 TFR (normalized)');

% Group 1
plot_spect_group1 = squeeze(nanmean(nanmean( ...
    grandavg_group1.powspctrm(:,:,7:179,36:171),2),1)); % freq × time
figure('Name','Group 1 TFR (MS\_F)');
imagesc(grandavg_group1.time(36:171), grandavg_group1.freq(7:179), plot_spect_group1);
set(gca,'YDir','normal');
colormap(jet);
ylim([4 90]);
xlim([-0.15 1.2]);
caxis([-0.5 0.5]);
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Grand average TFR – MS\_F');

% Group 2
plot_spect_group2 = squeeze(nanmean(nanmean( ...
    grandavg_group2.powspctrm(:,:,7:179,36:171),2),1));
figure('Name','Group 2 TFR (MS\_NF)');
imagesc(grandavg_group1.time(36:171), grandavg_group1.freq(7:179), plot_spect_group2);
set(gca,'YDir','normal');
colormap(jet);
ylim([4 90]);
xlim([-0.15 1.2]);
caxis([-0.5 0.5]);
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Grand average TFR – MS\_NF');

% Group 3
plot_spect_group3 = squeeze(nanmean(nanmean( ...
    grandavg_group3.powspctrm(:,:,7:179,36:171),2),1));
figure('Name','Group 3 TFR (HC)');
imagesc(grandavg_group1.time(36:171), grandavg_group1.freq(7:179), plot_spect_group3);
set(gca,'YDir','normal');
colormap(jet);
ylim([4 90]);
xlim([-0.15 1.2]);
caxis([-0.5 0.5]);
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Grand average TFR – HC');

