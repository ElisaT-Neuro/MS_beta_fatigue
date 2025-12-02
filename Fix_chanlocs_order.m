%% FIXCHANLOCS_TAPPING_MS
%
% Channel order correction for the MS tapping EEG dataset (FieldTrip format).
%
% Background:
%   In a subset of participants, a cabling error between the EEG cap and
%   the amplifier resulted in an incorrect order of channels at acquisition.
%   For these subjects, the data rows (channels) do not match the intended
%   61-channel 10–20 montage.
%
%   This script:
%     1) Identifies the subjects and conditions with the cabling error
%        (oldchanlocs list for each condition).
%     2) For those subjects:
%          - Reorders data.trial{t}(chan,:) so that the channels match the
%            canonical 61-channel montage.
%          - Reorders data.elec.elecpos, data.elec.pnt, and
%            data.elec.label consistently.
%          - Rewrites data.label to the same canonical order.
%     3) For subjects without cabling problems:
%          - Keeps the original data but ensures that labels and
%            electrode labels follow the canonical order.
%
%   The corrected FieldTrip structures are saved as:
%       Sbj<subject>_<cond>_fixchan.mat
%
%   These *_fixchan.mat files are then used for all subsequent analyses
%   (for example the TFR and source analyses) so that channel locations
%   and labels are consistent across all participants.
%
% Requirements:
%   - MATLAB
%   - FieldTrip-compatible data files:
%       SbjX_TAPPING_filt_notch_extrej_chanrej_ICA_IC_clean_ep_avref.mat
%   - labels_uppercase.mat: provides the canonical labels array "label"
%
% Script written by: Elisa Tatti
% Please do not reuse or adapt without permission.

clear; clc;

%% Base folders ----------------------------------------------------------

baseFolder = '/Users/leonardo/Dropbox (City College)/Tapping_EEG_MS_Siena/TAPPING_EEG/';
labelsFile = '/Users/leonardo/Dropbox (City College)/Tapping_EEG_MS_Siena/chanlocs/labels_uppercase.mat';

if ~exist(labelsFile, 'file')
    error('Labels file not found: %s', labelsFile);
end
load(labelsFile, 'label');    % canonical 61-channel label list

%% Group and subject definitions -----------------------------------------
% For each condition:
%   cond 1: MS_F (fatigued)
%   cond 2: MS_NF (non-fatigued)
%   cond 3: HC   (controls)
%
% "subjects" includes all subjects with usable EEG in that condition.
% "oldchanlocs" lists only those with the cabling error that require
% full reordering of channels and electrode positions.

for c = 1:3

    switch c
        case 1                                  % MS_F
            cond        = '1';
            subjects    = [1:9 11:19];          % 18 subjects; Sbj10 no EEG
            oldchanlocs = [1:9 11];             % subjects with cabling error

        case 2                                  % MS_NF
            cond        = '2';
            subjects    = [1 3:13 16:22];       % 20 subjects; Sbj15 no EEG; Sbj2 excluded
            oldchanlocs = [1:13];               % subjects with cabling error

        case 3                                  % HC
            cond        = '3';
            subjects    = [1:7 9:18];           % 17 subjects; Sbj8 wrong sampling rate
            oldchanlocs = [1 2 3:7];            % subjects with cabling error
    end

    folder_field = fullfile(baseFolder, cond, 'Preprocessing_avgref_ALE', 'Fieldtrip');

    fprintf('\n=== Condition %s ===\n', cond);
    fprintf('Subjects: %s\n', mat2str(subjects));
    fprintf('Subjects with cabling error (oldchanlocs): %s\n', mat2str(oldchanlocs));

    %% Loop over subjects ------------------------------------------------

    for s = subjects

        inFile = fullfile(folder_field, ...
            sprintf('Sbj%d_TAPPING_filt_notch_extrej_chanrej_ICA_IC_clean_ep_avref.mat', s));
        if ~exist(inFile, 'file')
            warning('Input file not found: %s (skipping Sbj%d, cond %s)', inFile, s, cond);
            continue;
        end

        fprintf('\nSubject %d, cond %s: loading %s\n', s, cond, inFile);
        load(inFile, 'data');   % FieldTrip data structure

        % Initialize data_new as a copy
        data_new = data;

        %% Reorder channels only if subject had cabling error -------------
        if ismember(s, oldchanlocs)

            fprintf('  -> Subject in oldchanlocs: applying channel reordering.\n');

            % ---------- TRIAL DATA REORDERING ----------
            % Each trial is a matrix [nChan x nTime].
            % We rewrite the rows in the canonical order.

            for t = 1:numel(data.trial)

                count = 1;

                % 1 Fp1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(1,:);
                count = count+1;

                % 2 Fpz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(58,:);
                count = count+1;

                % 3 Fp2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(2,:);
                count = count+1;

                % 4 AF3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(38,:);
                count = count+1;

                % 5 AF4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(39,:);
                count = count+1;

                % 6 F7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(11,:);
                count = count+1;

                % 7 F5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(46,:);
                count = count+1;

                % 8 F3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(3,:);
                count = count+1;

                % 9 F1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(32,:);
                count = count+1;

                % 10 Fz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(17,:);
                count = count+1;

                % 11 F2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(33,:);
                count = count+1;

                % 12 F4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(4,:);
                count = count+1;

                % 13 F6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(47,:);
                count = count+1;

                % 14 F8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(12,:);
                count = count+1;

                % 15 FT7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(52,:);
                count = count+1;

                % 16 FC5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(25,:);
                count = count+1;

                % 17 FC3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(40,:);
                count = count+1;

                % 18 FC1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(21,:);
                count = count+1;

                % 19 FCz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(20,:);
                count = count+1;

                % 20 FC2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(22,:);
                count = count+1;

                % 21 FC4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(41,:);
                count = count+1;

                % 22 FC6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(26,:);
                count = count+1;

                % 23 FT8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(53,:);
                count = count+1;

                % 24 T7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(13,:);
                count = count+1;

                % 25 C5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(48,:);
                count = count+1;

                % 26 C3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(5,:);
                count = count+1;

                % 27 C1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(34,:);
                count = count+1;

                % 28 Cz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(18,:);
                count = count+1;

                % 29 C2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(35,:);
                count = count+1;

                % 30 C4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(6,:);
                count = count+1;

                % 31 C6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(49,:);
                count = count+1;

                % 32 T8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(14,:);
                count = count+1;

                % 33 TP7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(54,:);
                count = count+1;

                % 34 CP5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(27,:);
                count = count+1;

                % 35 CP3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(42,:);
                count = count+1;

                % 36 CP1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(23,:);
                count = count+1;

                % 37 CPz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(59,:);
                count = count+1;

                % 38 CP2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(24,:);
                count = count+1;

                % 39 CP4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(43,:);
                count = count+1;

                % 40 CP6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(28,:);
                count = count+1;

                % 41 TP8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(55,:);
                count = count+1;

                % 42 P7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(15,:);
                count = count+1;

                % 43 P5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(50,:);
                count = count+1;

                % 44 P3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(7,:);
                count = count+1;

                % 45 P1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(36,:);
                count = count+1;

                % 46 Pz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(19,:);
                count = count+1;

                % 47 P2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(37,:);
                count = count+1;

                % 48 P4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(8,:);
                count = count+1;

                % 49 P6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(51,:);
                count = count+1;

                % 50 P8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(16,:);
                count = count+1;

                % 51 PO7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(56,:);
                count = count+1;

                % 52 PO3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(44,:);
                count = count+1;

                % 53 POz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(60,:);
                count = count+1;

                % 54 PO4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(45,:);
                count = count+1;

                % 55 PO8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(57,:);
                count = count+1;

                % 56 O1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(9,:);
                count = count+1;

                % 57 Oz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(61,:);
                count = count+1;

                % 58 O2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(10,:);
                count = count+1;

                % 59 TP9
                data_new.trial{1,t}(count,:) = data.trial{1,t}(30,:);
                count = count+1;

                % 60 TP10
                data_new.trial{1,t}(count,:) = data.trial{1,t}(31,:);
                count = count+1;

                % 61 AFz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(29,:);
                % count = count+1;  % not needed anymore
            end

            %% Electrode positions: elecpos --------------------------------

            count = 1;

            % Same mapping as above, but for data.elec.elecpos

            data_new.elec.elecpos(count,:) = data.elec.elecpos(1,:);   % Fp1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(58,:);  % Fpz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(2,:);   % Fp2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(38,:);  % AF3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(39,:);  % AF4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(11,:);  % F7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(46,:);  % F5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(3,:);   % F3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(32,:);  % F1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(17,:);  % Fz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(33,:);  % F2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(4,:);   % F4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(47,:);  % F6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(12,:);  % F8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(52,:);  % FT7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(25,:);  % FC5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(40,:);  % FC3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(21,:);  % FC1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(20,:);  % FCz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(22,:);  % FC2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(41,:);  % FC4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(26,:);  % FC6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(53,:);  % FT8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(13,:);  % T7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(48,:);  % C5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(5,:);   % C3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(34,:);  % C1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(18,:);  % Cz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(35,:);  % C2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(6,:);   % C4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(49,:);  % C6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(14,:);  % T8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(54,:);  % TP7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(27,:);  % CP5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(42,:);  % CP3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(23,:);  % CP1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(59,:);  % CPz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(24,:);  % CP2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(43,:);  % CP4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(28,:);  % CP6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(55,:);  % TP8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(15,:);  % P7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(50,:);  % P5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(7,:);   % P3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(36,:);  % P1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(19,:);  % Pz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(37,:);  % P2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(8,:);   % P4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(51,:);  % P6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(16,:);  % P8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(56,:);  % PO7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(44,:);  % PO3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(60,:);  % POz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(45,:);  % PO4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(57,:);  % PO8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(9,:);   % O1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(61,:);  % Oz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(10,:);  % O2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(30,:);  % TP9
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(31,:);  % TP10
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(29,:);  % AFz

            %% pnt (electrode positions in head coordinates) ----------------

            count = 1;
            data_new.elec.pnt(count,:) = data.elec.pnt(1,:);   % Fp1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(58,:);  % Fpz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(2,:);   % Fp2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(38,:);  % AF3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(39,:);  % AF4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(11,:);  % F7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(46,:);  % F5
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(3,:);   % F3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(32,:);  % F1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(17,:);  % Fz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(33,:);  % F2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(4,:);   % F4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(47,:);  % F6
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(12,:);  % F8
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(52,:);  % FT7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(25,:);  % FC5
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(40,:);  % FC3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(21,:);  % FC1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(20,:);  % FCz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(22,:);  % FC2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(41,:);  % FC4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(26,:);  % FC6
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(53,:);  % FT8
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(13,:);  % T7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(48,:);  % C5
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(5,:);   % C3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(34,:);  % C1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(18,:);  % Cz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(35,:);  % C2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(6,:);   % C4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(49,:);  % C6
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(14,:);  % T8
            count = count+1;
            data_new.elec._
%% FIXCHANLOCS_TAPPING_MS
%
% Channel order correction for the MS tapping EEG dataset (FieldTrip format).
%
% Background:
%   In a subset of participants, a cabling error between the EEG cap and
%   the amplifier resulted in an incorrect order of channels at acquisition.
%   For these subjects, the data rows (channels) do not match the intended
%   61-channel 10–20 montage.
%
%   This script:
%     1) Identifies the subjects and conditions with the cabling error
%        (oldchanlocs list for each condition).
%     2) For those subjects:
%          - Reorders data.trial{t}(chan,:) so that the channels match the
%            canonical 61-channel montage.
%          - Reorders data.elec.elecpos, data.elec.pnt, and
%            data.elec.label consistently.
%          - Rewrites data.label to the same canonical order.
%     3) For subjects without cabling problems:
%          - Keeps the original data but ensures that labels and
%            electrode labels follow the canonical order.
%
%   The corrected FieldTrip structures are saved as:
%       Sbj<subject>_<cond>_fixchan.mat
%
%   These *_fixchan.mat files are then used for all subsequent analyses
%   (for example the TFR and source analyses) so that channel locations
%   and labels are consistent across all participants.
%
% Requirements:
%   - MATLAB
%   - FieldTrip-compatible data files:
%       SbjX_TAPPING_filt_notch_extrej_chanrej_ICA_IC_clean_ep_avref.mat
%   - labels_uppercase.mat: provides the canonical labels array "label"
%
% Script written by: Elisa Tatti
% Please do not reuse or adapt without permission.

clear; clc;

%% Base folders ----------------------------------------------------------

baseFolder = '/Users/leonardo/Dropbox (City College)/Tapping_EEG_MS_Siena/TAPPING_EEG/';
labelsFile = '/Users/leonardo/Dropbox (City College)/Tapping_EEG_MS_Siena/chanlocs/labels_uppercase.mat';

if ~exist(labelsFile, 'file')
    error('Labels file not found: %s', labelsFile);
end
load(labelsFile, 'label');    % canonical 61-channel label list

%% Group and subject definitions -----------------------------------------
% For each condition:
%   cond 1: MS_F (fatigued)
%   cond 2: MS_NF (non-fatigued)
%   cond 3: HC   (controls)
%
% "subjects" includes all subjects with usable EEG in that condition.
% "oldchanlocs" lists only those with the cabling error that require
% full reordering of channels and electrode positions.

for c = 1:3

    switch c
        case 1                                  % MS_F
            cond        = '1';
            subjects    = [1:9 11:19];          % 18 subjects; Sbj10 no EEG
            oldchanlocs = [1:9 11];             % subjects with cabling error

        case 2                                  % MS_NF
            cond        = '2';
            subjects    = [1 3:13 16:22];       % 20 subjects; Sbj15 no EEG; Sbj2 excluded
            oldchanlocs = [1:13];               % subjects with cabling error

        case 3                                  % HC
            cond        = '3';
            subjects    = [1:7 9:18];           % 17 subjects; Sbj8 wrong sampling rate
            oldchanlocs = [1 2 3:7];            % subjects with cabling error
    end

    folder_field = fullfile(baseFolder, cond, 'Preprocessing_avgref_ALE', 'Fieldtrip');

    fprintf('\n=== Condition %s ===\n', cond);
    fprintf('Subjects: %s\n', mat2str(subjects));
    fprintf('Subjects with cabling error (oldchanlocs): %s\n', mat2str(oldchanlocs));

    %% Loop over subjects ------------------------------------------------

    for s = subjects

        inFile = fullfile(folder_field, ...
            sprintf('Sbj%d_TAPPING_filt_notch_extrej_chanrej_ICA_IC_clean_ep_avref.mat', s));
        if ~exist(inFile, 'file')
            warning('Input file not found: %s (skipping Sbj%d, cond %s)', inFile, s, cond);
            continue;
        end

        fprintf('\nSubject %d, cond %s: loading %s\n', s, cond, inFile);
        load(inFile, 'data');   % FieldTrip data structure

        % Initialize data_new as a copy
        data_new = data;

        %% Reorder channels only if subject had cabling error -------------
        if ismember(s, oldchanlocs)

            fprintf('  -> Subject in oldchanlocs: applying channel reordering.\n');

            % ---------- TRIAL DATA REORDERING ----------
            % Each trial is a matrix [nChan x nTime].
            % We rewrite the rows in the canonical order.

            for t = 1:numel(data.trial)

                count = 1;

                % 1 Fp1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(1,:);
                count = count+1;

                % 2 Fpz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(58,:);
                count = count+1;

                % 3 Fp2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(2,:);
                count = count+1;

                % 4 AF3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(38,:);
                count = count+1;

                % 5 AF4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(39,:);
                count = count+1;

                % 6 F7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(11,:);
                count = count+1;

                % 7 F5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(46,:);
                count = count+1;

                % 8 F3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(3,:);
                count = count+1;

                % 9 F1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(32,:);
                count = count+1;

                % 10 Fz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(17,:);
                count = count+1;

                % 11 F2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(33,:);
                count = count+1;

                % 12 F4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(4,:);
                count = count+1;

                % 13 F6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(47,:);
                count = count+1;

                % 14 F8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(12,:);
                count = count+1;

                % 15 FT7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(52,:);
                count = count+1;

                % 16 FC5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(25,:);
                count = count+1;

                % 17 FC3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(40,:);
                count = count+1;

                % 18 FC1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(21,:);
                count = count+1;

                % 19 FCz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(20,:);
                count = count+1;

                % 20 FC2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(22,:);
                count = count+1;

                % 21 FC4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(41,:);
                count = count+1;

                % 22 FC6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(26,:);
                count = count+1;

                % 23 FT8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(53,:);
                count = count+1;

                % 24 T7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(13,:);
                count = count+1;

                % 25 C5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(48,:);
                count = count+1;

                % 26 C3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(5,:);
                count = count+1;

                % 27 C1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(34,:);
                count = count+1;

                % 28 Cz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(18,:);
                count = count+1;

                % 29 C2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(35,:);
                count = count+1;

                % 30 C4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(6,:);
                count = count+1;

                % 31 C6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(49,:);
                count = count+1;

                % 32 T8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(14,:);
                count = count+1;

                % 33 TP7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(54,:);
                count = count+1;

                % 34 CP5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(27,:);
                count = count+1;

                % 35 CP3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(42,:);
                count = count+1;

                % 36 CP1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(23,:);
                count = count+1;

                % 37 CPz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(59,:);
                count = count+1;

                % 38 CP2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(24,:);
                count = count+1;

                % 39 CP4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(43,:);
                count = count+1;

                % 40 CP6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(28,:);
                count = count+1;

                % 41 TP8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(55,:);
                count = count+1;

                % 42 P7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(15,:);
                count = count+1;

                % 43 P5
                data_new.trial{1,t}(count,:) = data.trial{1,t}(50,:);
                count = count+1;

                % 44 P3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(7,:);
                count = count+1;

                % 45 P1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(36,:);
                count = count+1;

                % 46 Pz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(19,:);
                count = count+1;

                % 47 P2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(37,:);
                count = count+1;

                % 48 P4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(8,:);
                count = count+1;

                % 49 P6
                data_new.trial{1,t}(count,:) = data.trial{1,t}(51,:);
                count = count+1;

                % 50 P8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(16,:);
                count = count+1;

                % 51 PO7
                data_new.trial{1,t}(count,:) = data.trial{1,t}(56,:);
                count = count+1;

                % 52 PO3
                data_new.trial{1,t}(count,:) = data.trial{1,t}(44,:);
                count = count+1;

                % 53 POz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(60,:);
                count = count+1;

                % 54 PO4
                data_new.trial{1,t}(count,:) = data.trial{1,t}(45,:);
                count = count+1;

                % 55 PO8
                data_new.trial{1,t}(count,:) = data.trial{1,t}(57,:);
                count = count+1;

                % 56 O1
                data_new.trial{1,t}(count,:) = data.trial{1,t}(9,:);
                count = count+1;

                % 57 Oz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(61,:);
                count = count+1;

                % 58 O2
                data_new.trial{1,t}(count,:) = data.trial{1,t}(10,:);
                count = count+1;

                % 59 TP9
                data_new.trial{1,t}(count,:) = data.trial{1,t}(30,:);
                count = count+1;

                % 60 TP10
                data_new.trial{1,t}(count,:) = data.trial{1,t}(31,:);
                count = count+1;

                % 61 AFz
                data_new.trial{1,t}(count,:) = data.trial{1,t}(29,:);
                % count = count+1;  % not needed anymore
            end

            %% Electrode positions: elecpos --------------------------------

            count = 1;

            % Same mapping as above, but for data.elec.elecpos

            data_new.elec.elecpos(count,:) = data.elec.elecpos(1,:);   % Fp1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(58,:);  % Fpz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(2,:);   % Fp2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(38,:);  % AF3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(39,:);  % AF4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(11,:);  % F7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(46,:);  % F5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(3,:);   % F3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(32,:);  % F1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(17,:);  % Fz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(33,:);  % F2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(4,:);   % F4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(47,:);  % F6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(12,:);  % F8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(52,:);  % FT7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(25,:);  % FC5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(40,:);  % FC3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(21,:);  % FC1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(20,:);  % FCz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(22,:);  % FC2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(41,:);  % FC4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(26,:);  % FC6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(53,:);  % FT8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(13,:);  % T7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(48,:);  % C5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(5,:);   % C3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(34,:);  % C1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(18,:);  % Cz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(35,:);  % C2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(6,:);   % C4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(49,:);  % C6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(14,:);  % T8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(54,:);  % TP7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(27,:);  % CP5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(42,:);  % CP3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(23,:);  % CP1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(59,:);  % CPz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(24,:);  % CP2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(43,:);  % CP4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(28,:);  % CP6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(55,:);  % TP8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(15,:);  % P7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(50,:);  % P5
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(7,:);   % P3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(36,:);  % P1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(19,:);  % Pz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(37,:);  % P2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(8,:);   % P4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(51,:);  % P6
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(16,:);  % P8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(56,:);  % PO7
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(44,:);  % PO3
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(60,:);  % POz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(45,:);  % PO4
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(57,:);  % PO8
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(9,:);   % O1
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(61,:);  % Oz
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(10,:);  % O2
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(30,:);  % TP9
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(31,:);  % TP10
            count = count+1;
            data_new.elec.elecpos(count,:) = data.elec.elecpos(29,:);  % AFz

            %% pnt (electrode positions in head coordinates) ----------------

            count = 1;
            data_new.elec.pnt(count,:) = data.elec.pnt(1,:);   % Fp1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(58,:);  % Fpz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(2,:);   % Fp2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(38,:);  % AF3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(39,:);  % AF4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(11,:);  % F7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(46,:);  % F5
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(3,:);   % F3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(32,:);  % F1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(17,:);  % Fz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(33,:);  % F2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(4,:);   % F4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(47,:);  % F6
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(12,:);  % F8
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(52,:);  % FT7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(25,:);  % FC5
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(40,:);  % FC3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(21,:);  % FC1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(20,:);  % FCz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(22,:);  % FC2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(41,:);  % FC4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(26,:);  % FC6
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(53,:);  % FT8
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(13,:);  % T7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(48,:);  % C5
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(5,:);   % C3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(34,:);  % C1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(18,:);  % Cz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(35,:);  % C2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(6,:);   % C4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(49,:);  % C6
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(14,:);  % T8
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(54,:);  % TP7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(27,:);  % CP5
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(42,:);  % CP3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(23,:);  % CP1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(59,:);  % CPz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(24,:);  % CP2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(43,:);  % CP4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(28,:);  % CP6
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(55,:);  % TP8
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(15,:);  % P7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(50,:);  % P5
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(7,:);   % P3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(36,:);  % P1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(19,:);  % Pz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(37,:);  % P2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(8,:);   % P4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(51,:);  % P6
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(16,:);  % P8
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(56,:);  % PO7
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(44,:);  % PO3
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(60,:);  % POz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(45,:);  % PO4
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(57,:);  % PO8
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(9,:);   % O1
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(61,:);  % Oz
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(10,:);  % O2
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(30,:);  % TP9
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(31,:);  % TP10
            count = count+1;
            data_new.elec.pnt(count,:) = data.elec.pnt(29,:);  % AFz
        else
            fprintf('  -> Subject not in oldchanlocs: channel order kept, only labels enforced.\n');
        end

        %% Canonical labels for all subjects ---------------------------------
        % For both corrected and non corrected subjects, enforce the canonical
        % label order that is stored in labels_uppercase.mat

        data_new.label       = label;
        data_new.elec.label  = label;

        % You can optionally check the layout here:
        % cfg = [];
        % layout = ft_prepare_layout(cfg, data_new);
        % figure; ft_plot_layout(layout);

        %% Save corrected data ----------------------------------------------

        outFile = fullfile(folder_field, sprintf('Sbj%d_%s_fixchan.mat', s, cond));
        data    = data_new; %#ok<NASGU>
        save(outFile, 'data', '-v7.3');
        fprintf('  -> Saved corrected data to %s\n', outFile);

        clear data data_new;
    end
end
