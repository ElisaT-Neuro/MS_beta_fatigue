function extract_beta_modulation()
% EXTRACT_BETA_MODULATION
%
% Computes movement-related beta ERD/ERS and modulation depth from single-trial
% Laplacian-filtered EEG using a Morlet wavelet time–frequency analysis.
% Outputs single-trial and trial-averaged measures for left, right and frontal
% ROIs as .csv tables.
%
% This code was used in:
%   Tatti et al., "Insights into Central Fatigue in Multiple Sclerosis:
%   Movement-Related Beta Oscillatory Activity and its Association with
%   Brain Neurophysiology and Structure" (Brain Communications).
%
% REQUIREMENTS
%   - MATLAB R20xx
%   - FieldTrip on your MATLAB path (https://www.fieldtriptoolbox.org)
%   - Preprocessed EEG data saved as FieldTrip structures:
%       data.elec  : electrode positions
%       data.trial : trials x channels x time
%       data.time  : 1 x time vector
%
%   - Excel files with ROIs for each group, with sheets 'left', 'right', 'frontal'.
%
% NOTE
%   - This script is written as a single function for transparency.
%   - To adapt it to your own project, change only the paths and subject lists
%     in the USER SETTINGS section below.

%% ====================== USER SETTINGS ==================================

% Project root (edit this to point to your project folder)
projectRoot = '/path/to/TAPPING_EEG_MS_Siena/TAPPING_EEG';  % <-- CHANGE THIS

% Folder containing preprocessed FieldTrip data:
% expected structure: <projectRoot>/<cond>/Preprocessing_avgref_ALE/Fieldtrip/Sbj<s>_<cond>_fixchan.mat
dataRoot = projectRoot;

% Output folders (they will be created if they do not exist)
outputRoot        = fullfile(projectRoot, 'Analyses_beta_modulation_NewTFR');
tfrFolder         = fullfile(outputRoot, 'TFRs');
singleTrialFolder = fullfile(outputRoot, 'Tables_betamod_singletrial');
timeFolder        = fullfile(outputRoot, 'Tables_betamod_time');
avgPowFolder      = fullfile(outputRoot, 'Tables_avg_csv', 'Avg_trials_betamodulation');

if ~exist(outputRoot, 'dir');        mkdir(outputRoot);        end
if ~exist(tfrFolder, 'dir');         mkdir(tfrFolder);         end
if ~exist(singleTrialFolder, 'dir'); mkdir(singleTrialFolder); end
if ~exist(timeFolder, 'dir');        mkdir(timeFolder);        end
if ~exist(avgPowFolder, 'dir');      mkdir(avgPowFolder);      end

% Group/condition definitions
% cond codes must match the subfolder names under projectRoot (e.g. '1','2','3').
groups(1).name        = 'MS_F';
groups(1).cond        = '1';
groups(1).subjects    = [1 2 3 4:9 11:19];  % example list from original script
groups(1).roiFilePath = fullfile(projectRoot, 'Analyses_beta_modulation', 'ROIs_betamod_Group1.xlsx');

groups(2).name        = 'MS_NF';
groups(2).cond        = '2';
groups(2).subjects    = [1 2 3:14 16:22];   % Sbj15 missing; Sbj2 very noisy
groups(2).roiFilePath = fullfile(projectRoot, 'Analyses_beta_modulation', 'ROIs_betamod_Group2.xlsx');

groups(3).name        = 'HC';
groups(3).cond        = '3';
groups(3).subjects    = [1 2 3:7 9:18];     % Sbj8 deleted (wrong sampling rate)
groups(3).roiFilePath = fullfile(projectRoot, 'Analyses_beta_modulation', 'ROIs_betamod_Group3.xlsx');

% Which groups to run (indices into 'groups')
groupsToRun = [2];  % e.g. [1 2 3] to run all; here 2 = MS_NF as in your example

% Preallocation sizes (adapt if needed)
maxTrials   = 310;  % max number of trials per subject
maxSubjects = 25;   % max number of subjects per group

%% =================== LOOP OVER GROUPS / CONDITIONS =====================

for g = groupsToRun

    groupName  = groups(g).name;
    cond       = groups(g).cond;
    subjects   = groups(g).subjects;
    roiFile    = groups(g).roiFilePath;

    fprintf('Processing group %s (cond %s)\n', groupName, cond);

    % Load ROIs for this group
    roi_right   = xlsread(roiFile, 'right');
    roi_left    = xlsread(roiFile, 'left');
    roi_frontal = xlsread(roiFile, 'frontal');

    % Preallocate group-level matrices
    Left_Brd    = nan(maxTrials, maxSubjects);
    Left_Brs    = nan(maxTrials, maxSubjects);
    Left_Bp2p   = nan(maxTrials, maxSubjects);

    Right_Brd   = nan(maxTrials, maxSubjects);
    Right_Brs   = nan(maxTrials, maxSubjects);
    Right_Bp2p  = nan(maxTrials, maxSubjects);

    Frontal_Brd  = nan(maxTrials, maxSubjects);
    Frontal_Brs  = nan(maxTrials, maxSubjects);
    Frontal_Bp2p = nan(maxTrials, maxSubjects);

    Brd_left_time_avg   = nan(maxSubjects, 1);
    Brs_left_time_avg   = nan(maxSubjects, 1);
    Brd_right_time_avg  = nan(maxSubjects, 1);
    Brs_right_time_avg  = nan(maxSubjects, 1);
    Brd_front_time_avg  = nan(maxSubjects, 1);
    Brs_front_time_avg  = nan(maxSubjects, 1);

    Brd_left_pow_avg    = nan(maxSubjects, 1);
    Brs_left_pow_avg    = nan(maxSubjects, 1);
    Bp2p_left_pow_avg   = nan(maxSubjects, 1);

    Brd_right_pow_avg   = nan(maxSubjects, 1);
    Brs_right_pow_avg   = nan(maxSubjects, 1);
    Bp2p_right_pow_avg  = nan(maxSubjects, 1);

    Brd_front_pow_avg   = nan(maxSubjects, 1);
    Brs_front_pow_avg   = nan(maxSubjects, 1);
    Bp2p_front_pow_avg  = nan(maxSubjects, 1);

    % Time-at-peak for each trial and subject
    Brd_left_time   = nan(maxTrials, maxSubjects);
    Brs_left_time   = nan(maxTrials, maxSubjects);
    Brd_right_time  = nan(maxTrials, maxSubjects);
    Brs_right_time  = nan(maxTrials, maxSubjects);
    Brd_front_time  = nan(maxTrials, maxSubjects);
    Brs_front_time  = nan(maxTrials, maxSubjects);

    counts = 1;  % subject counter within group for indexing columns

    %% ============== LOOP OVER SUBJECTS =================================

    for s = subjects

        fprintf('  Subject %d\n', s);

        % Load preprocessed data (average reference, laplacian will be applied here)
        dataFile = fullfile(dataRoot, cond, 'Preprocessing_avgref_ALE', 'Fieldtrip', ...
                            sprintf('Sbj%d_%s_fixchan.mat', s, cond));
        if ~exist(dataFile, 'file')
            warning('File not found: %s. Skipping subject.', dataFile);
            continue;
        end
        load(dataFile, 'data');  % loads variable 'data' (FieldTrip structure)

        %% Prepare layout and apply Laplacian filter
        cfg            = [];
        cfg.elec       = data.elec;
        cfg.rotate     = [];              % default projection
        cfg.projection = 'orthographic';
        layout         = ft_prepare_layout(cfg, data);

        cfg            = [];
        cfg.method     = 'laplacian';
        cfg.layout     = layout;
        cfg.feedback   = 'no';
        data_lap       = ft_preprocessing(cfg, data);

        %% Time–frequency analysis using Morlet wavelets

        f_o_i   = 1:0.05:90;             % frequencies of interest (Hz)
        t_o_i   = -0.5:0.01:1.5;         % time of interest (s)

        cfg            = [];
        cfg.method     = 'wavelet';
        cfg.output     = 'pow';
        cfg.channel    = 'all';
        cfg.trials     = 'all';
        cfg.keeptrials = 'yes';
        cfg.keeptapers = 'no';
        cfg.foi        = f_o_i;
        cfg.toi        = t_o_i;
        cfg.pad        = 4;              % zero-padding to increase frequency resolution

        % Frequency-dependent wavelet width (3 to 10 cycles)
        minCycles      = 3;
        maxCycles      = 10;
        cfg.width      = linspace(minCycles, maxCycles, numel(f_o_i));

        freq = ft_freqanalysis(cfg, data_lap);
        save(fullfile(tfrFolder, sprintf('Sbj_%d_%s_TFR.mat', s, groupName)), 'freq');

        %% Normalize power (trial-wise, relative change from baseline)
        % freq.powspctrm: trials x channels x frequencies x time

        freq_norm = freq;      % keep structure, replace powspctrm
        dims      = size(freq_norm.powspctrm);
        trials    = dims(1);

        % Baseline window: here from time indices 1:201 (adapt if needed)
        % refblock: chan x freq
        refblock = squeeze(nanmean(nanmean(freq_norm.powspctrm(:,:,:,1:201), 1), 4));

        % Repeat baseline across time dimension
        nTime      = 201; % corresponding to 0–2 s in original script
        refblock_t = repmat(refblock, [1 1 nTime]);

        normtrials = nan(trials, dims(2), dims(3), nTime);
        for t = 1:trials
            normtrials(t,:,:,:) = (squeeze(freq_norm.powspctrm(t,:,:,:)) - refblock_t) ./ refblock_t;
        end

        freq_norm.powspctrm = normtrials;
        clear freq normtrials refblock refblock_t;

        %% Helper: beta-band average across frequencies of interest
        % original script uses indices 51:97 in freq dimension (beta band)
        betaFreqIdx = 51:97;

        %% LEFT ROI
        roi_l = roi_left(s, 3:9);  % channels for this subject
        if any(~isnan(roi_l))
            [Left_Brd, Left_Brs, Left_Bp2p, ...
             Brd_left_time, Brs_left_time, ...
             Brd_left_time_avg, Brs_left_time_avg, ...
             Brd_left_pow_avg, Brs_left_pow_avg, Bp2p_left_pow_avg] = ...
                 compute_roi_beta_modulation(freq_norm, roi_l, betaFreqIdx, ...
                                             Left_Brd, Left_Brs, Left_Bp2p, ...
                                             Brd_left_time, Brs_left_time, ...
                                             Brd_left_time_avg, Brs_left_time_avg, ...
                                             Brd_left_pow_avg, Brs_left_pow_avg, Bp2p_left_pow_avg, ...
                                             counts);
        else
            Brd_left_time_avg(counts,:) = NaN;
            Brs_left_time_avg(counts,:) = NaN;
        end

        %% RIGHT ROI
        roi_r = roi_right(s, 3:9);
        if any(~isnan(roi_r))
            [Right_Brd, Right_Brs, Right_Bp2p, ...
             Brd_right_time, Brs_right_time, ...
             Brd_right_time_avg, Brs_right_time_avg, ...
             Brd_right_pow_avg, Brs_right_pow_avg, Bp2p_right_pow_avg] = ...
                 compute_roi_beta_modulation(freq_norm, roi_r, betaFreqIdx, ...
                                             Right_Brd, Right_Brs, Right_Bp2p, ...
                                             Brd_right_time, Brs_right_time, ...
                                             Brd_right_time_avg, Brs_right_time_avg, ...
                                             Brd_right_pow_avg, Brs_right_pow_avg, Bp2p_right_pow_avg, ...
                                             counts);
        else
            Brd_right_time_avg(counts,:) = NaN;
            Brs_right_time_avg(counts,:) = NaN;
        end

        %% FRONTAL ROI
        roi_f = roi_frontal(s, 3:9);
        if any(~isnan(roi_f))
            [Frontal_Brd, Frontal_Brs, Frontal_Bp2p, ...
             Brd_front_time, Brs_front_time, ...
             Brd_front_time_avg, Brs_front_time_avg, ...
             Brd_front_pow_avg, Brs_front_pow_avg, Bp2p_front_pow_avg] = ...
                 compute_roi_beta_modulation(freq_norm, roi_f, betaFreqIdx, ...
                                             Frontal_Brd, Frontal_Brs, Frontal_Bp2p, ...
                                             Brd_front_time, Brs_front_time, ...
                                             Brd_front_time_avg, Brs_front_time_avg, ...
                                             Brd_front_pow_avg, Brs_front_pow_avg, Bp2p_front_pow_avg, ...
                                             counts);
        else
            Brd_front_time_avg(counts,:) = NaN;
            Brs_front_time_avg(counts,:) = NaN;
        end

        counts = counts + 1;

        clear data data_lap freq_norm roi_l roi_r roi_f;
    end % subjects

    %% ================= SAVE TABLES FOR THIS GROUP =====================

    condStr = cond;  % for filenames

    % Single-trial amplitudes (ERD, ERS, peak-to-peak) by ROI

    Left_Brd_table   = array2table(Left_Brd);
    Left_Brs_table   = array2table(Left_Brs);
    Left_Bp2p_table  = array2table(Left_Bp2p);

    Right_Brd_table  = array2table(Right_Brd);
    Right_Brs_table  = array2table(Right_Brs);
    Right_Bp2p_table = array2table(Right_Bp2p);

    Front_Brd_table  = array2table(Frontal_Brd);
    Front_Brs_table  = array2table(Frontal_Brs);
    Front_Bp2p_table = array2table(Frontal_Bp2p);

    % Time of ERD/ERS peaks (all trials and trial averages)

    Brd_left_time_table   = array2table(Brd_left_time);
    Brs_left_time_table   = array2table(Brs_left_time);
    Brd_right_time_table  = array2table(Brd_right_time);
    Brs_right_time_table  = array2table(Brs_right_time);
    Brd_front_time_table  = array2table(Brd_front_time);
    Brs_front_time_table  = array2table(Brs_front_time);

    Brd_left_time_avg_table   = array2table(Brd_left_time_avg);
    Brs_left_time_avg_table   = array2table(Brs_left_time_avg);
    Brd_right_time_avg_table  = array2table(Brd_right_time_avg);
    Brs_right_time_avg_table  = array2table(Brs_right_time_avg);
    Brd_front_time_avg_table  = array2table(Brd_front_time_avg);
    Brs_front_time_avg_table  = array2table(Brs_front_time_avg);

    % Save time (all trials)
    writetable(Brd_left_time_table,  fullfile(timeFolder, sprintf('Table_MS_%s_Left_Brd_time_all_trials.csv',   condStr)));
    writetable(Brs_left_time_table,  fullfile(timeFolder, sprintf('Table_MS_%s_Left_Brs_time_all_trials.csv',   condStr)));
    writetable(Brd_front_time_table, fullfile(timeFolder, sprintf('Table_MS_%s_Front_Brd_time_all_trials.csv',  condStr)));
    writetable(Brs_front_time_table, fullfile(timeFolder, sprintf('Table_MS_%s_Front_Brs_time_all_trials.csv',  condStr)));
    writetable(Brd_right_time_table, fullfile(timeFolder, sprintf('Table_MS_%s_Right_Brd_time_all_trials.csv',  condStr)));
    writetable(Brs_right_time_table, fullfile(timeFolder, sprintf('Table_MS_%s_Right_Brs_time_all_trials.csv',  condStr)));

    % Save time (average across trials)
    writetable(Brd_left_time_avg_table,  fullfile(timeFolder, sprintf('Table_MS_%s_Left_Brd_time_avg_trials.csv',   condStr)));
    writetable(Brs_left_time_avg_table,  fullfile(timeFolder, sprintf('Table_MS_%s_Left_Brs_time_avg_trials.csv',   condStr)));
    writetable(Brd_right_time_avg_table, fullfile(timeFolder, sprintf('Table_MS_%s_Right_Brd_time_avg_trials.csv',  condStr)));
    writetable(Brs_right_time_avg_table, fullfile(timeFolder, sprintf('Table_MS_%s_Right_Brs_time_avg_trials.csv',  condStr)));
    writetable(Brd_front_time_avg_table, fullfile(timeFolder, sprintf('Table_MS_%s_Front_Brd_time_avg_trials.csv',  condStr)));
    writetable(Brs_front_time_avg_table, fullfile(timeFolder, sprintf('Table_MS_%s_Front_Brs_time_avg_trials.csv',  condStr)));

    % Save single-trial amplitudes
    writetable(Left_Brs_table,   fullfile(singleTrialFolder, sprintf('Table_MS_%s_Left_Brs_RM_singletrial.csv',   condStr)));
    writetable(Left_Brd_table,   fullfile(singleTrialFolder, sprintf('Table_MS_%s_Left_Brd_RM_singletrial.csv',   condStr)));
    writetable(Left_Bp2p_table,  fullfile(singleTrialFolder, sprintf('Table_MS_%s_Left_Bp2p_RM_singletrial.csv',  condStr)));

    writetable(Right_Bp2p_table, fullfile(singleTrialFolder, sprintf('Table_MS_%s_Right_Bp2p_RM_singletrial.csv', condStr)));
    writetable(Right_Brs_table,  fullfile(singleTrialFolder, sprintf('Table_MS_%s_Right_Brs_RM_singletrial.csv',  condStr)));
    writetable(Right_Brd_table,  fullfile(singleTrialFolder, sprintf('Table_MS_%s_Right_Brd_RM_singletrial.csv',  condStr)));

    writetable(Front_Bp2p_table, fullfile(singleTrialFolder, sprintf('Table_MS_%s_Front_Bp2p_RM_singletrial.csv', condStr)));
    writetable(Front_Brs_table,  fullfile(singleTrialFolder, sprintf('Table_MS_%s_Front_Brs_RM_singletrial.csv',  condStr)));
    writetable(Front_Brd_table,  fullfile(singleTrialFolder, sprintf('Table_MS_%s_Front_Brd_RM_singletrial.csv',  condStr)));

    % Trial-averaged amplitudes (note: fixed the name swap bug here)
    Left_Bp2p_avg_table  = array2table(Bp2p_left_pow_avg);
    Left_Brd_avg_table   = array2table(Brd_left_pow_avg);
    Left_Brs_avg_table   = array2table(Brs_left_pow_avg);

    Right_Bp2p_avg_table = array2table(Bp2p_right_pow_avg);
    Right_Brd_avg_table  = array2table(Brd_right_pow_avg);
    Right_Brs_avg_table  = array2table(Brs_right_pow_avg);

    Front_Bp2p_avg_table = array2table(Bp2p_front_pow_avg);
    Front_Brd_avg_table  = array2table(Brd_front_pow_avg);
    Front_Brs_avg_table  = array2table(Brs_front_pow_avg);

    % Save averaged amplitudes
    writetable(Left_Bp2p_avg_table,  fullfile(avgPowFolder, sprintf('Table_MS_%s_Left_Bp2p_avg_trials.csv',  condStr)));
    writetable(Left_Brd_avg_table,   fullfile(avgPowFolder, sprintf('Table_MS_%s_Left_Brd_avg_trials.csv',   condStr)));
    writetable(Left_Brs_avg_table,   fullfile(avgPowFolder, sprintf('Table_MS_%s_Left_Brs_avg_trials.csv',   condStr)));

    writetable(Right_Bp2p_avg_table, fullfile(avgPowFolder, sprintf('Table_MS_%s_Right_Bp2p_avg_trials.csv', condStr)));
    writetable(Right_Brs_avg_table,  fullfile(avgPowFolder, sprintf('Table_MS_%s_Right_Brs_avg_trials.csv',  condStr)));
    writetable(Right_Brd_avg_table,  fullfile(avgPowFolder, sprintf('Table_MS_%s_Right_Brd_avg_trials.csv',  condStr)));

    writetable(Front_Bp2p_avg_table, fullfile(avgPowFolder, sprintf('Table_MS_%s_Front_Bp2p_avg_trials.csv', condStr)));
    writetable(Front_Brs_avg_table,  fullfile(avgPowFolder, sprintf('Table_MS_%s_Front_Brs_avg_trials.csv',  condStr)));
    writetable(Front_Brd_avg_table,  fullfile(avgPowFolder, sprintf('Table_MS_%s_Front_Brd_avg_trials.csv',  condStr)));

    fprintf('Finished group %s (cond %s)\n', groupName, condStr);
end

end % main function


function [ROI_Brd, ROI_Brs, ROI_Bp2p, ...
          Brd_time_mat, Brs_time_mat, ...
          Brd_time_avg, Brs_time_avg, ...
          Brd_pow_avg, Brs_pow_avg, Bp2p_pow_avg] = ...
          compute_roi_beta_modulation(freq_norm, roi_idx, betaFreqIdx, ...
                                      ROI_Brd, ROI_Brs, ROI_Bp2p, ...
                                      Brd_time_mat, Brs_time_mat, ...
                                      Brd_time_avg, Brs_time_avg, ...
                                      Brd_pow_avg, Brs_pow_avg, Bp2p_pow_avg, ...
                                      subjIndex)
% Helper function to compute ERD, ERS and peak-to-peak modulation for a given ROI.

trials = size(freq_norm.powspctrm, 1);

% Average beta power across frequencies in the beta band
Bnormtrials = squeeze(mean(freq_norm.powspctrm(:,:,betaFreqIdx,:), 3));   % trials x channels x time
Bnormtrials = squeeze(mean(Bnormtrials(:, roi_idx, :), 2));               % trials x time

% Find ERD (minimum) and ERS (maximum) in predefined time windows
% time indices based on original script:
% ERD: indices 61:121   (~0.1–0.7 s)
% ERS: indices 121:191  (~1–2 s)
[~, Brd_idx] = min(Bnormtrials(:,61:121), [], 2);
[~, Brs_idx] = max(Bnormtrials(:,121:191), [], 2);

BrdT = Brd_idx + 60;
BrsT = Brs_idx + 120;

% Extract time–frequency slices at ERD and ERS peak times
Brd = nan(trials, size(freq_norm.powspctrm, 2), size(freq_norm.powspctrm, 3));
Brs = nan(trials, size(freq_norm.powspctrm, 2), size(freq_norm.powspctrm, 3));

Brd_t = nan(trials, 1);
Brs_t = nan(trials, 1);

for t = 1:trials
    Brd(t,:,:) = squeeze(freq_norm.powspctrm(t,:,:,BrdT(t)));
    Brs(t,:,:) = squeeze(freq_norm.powspctrm(t,:,:,BrsT(t)));

    Brd_t(t)   = freq_norm.time(BrdT(t));
    Brs_t(t)   = freq_norm.time(BrsT(t));
end

Bp2p = Brs - Brd;

% Average within ROI and beta band
Brd_roi   = squeeze(nanmean(nanmean(Brd(:, roi_idx, betaFreqIdx), 3), 2));
Brs_roi   = squeeze(nanmean(nanmean(Brs(:, roi_idx, betaFreqIdx), 3), 2));
Bp2p_roi  = squeeze(nanmean(nanmean(Bp2p(:, roi_idx, betaFreqIdx), 3), 2));

% Fill group-level matrices for this subject (column = subjIndex)
ROI_Brd(1:trials, subjIndex)   = Brd_roi;
ROI_Brs(1:trials, subjIndex)   = Brs_roi;
ROI_Bp2p(1:trials, subjIndex)  = Bp2p_roi;

Brd_time_mat(1:trials, subjIndex) = Brd_t;
Brs_time_mat(1:trials, subjIndex) = Brs_t;

Brd_time_avg(subjIndex,:) = nanmean(Brd_t, 1);
Brs_time_avg(subjIndex,:) = nanmean(Brs_t, 1);

Brd_pow_avg(subjIndex,:)  = nanmean(Brd_roi, 1);
Brs_pow_avg(subjIndex,:)  = nanmean(Brs_roi, 1);
Bp2p_pow_avg(subjIndex,:) = nanmean(Bp2p_roi, 1);

end
