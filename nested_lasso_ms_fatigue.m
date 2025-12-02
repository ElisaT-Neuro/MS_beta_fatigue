function nested_lasso_ms_fatigue()
% NESTED_LASSO_MS_FATIGUE
%
% Nested elastic-net / LASSO logistic regression for MS fatigue classification.
%
%  - Outcome coding (IMPORTANT):
%       y = 1 -> fatigued (MS_F)
%       y = 0 -> non-fatigued (MS_NF)
%  - Theory-driven candidate predictors; NO outcome-based pre-screening
%  - Stratified Kouter x Kinner nested CV for unbiased performance
%  - Repeated nested CV for confidence intervals and selection stability
%  - Descriptive penalized (near-ridge) refit on features with selection
%    frequency >= 0.50
%  - Reports β, odds ratios (per 1 SD), 95% bootstrap CIs, and pseudo-R²
%
% This script was used in:
%   Tatti et al., "Insights into Central Fatigue in Multiple Sclerosis:
%   Movement-Related Beta Oscillatory Activity and its Association with
%   Brain Neurophysiology and Structure" (Brain Communications).
%
% DEPENDENCIES
%   - Statistics and Machine Learning Toolbox (for lassoglm, cvpartition, perfcurve)
%
% USAGE (typical):
%   1. Place your analysis CSV in:  data/lasso_with_caudateandSLF.csv
%   2. Edit the 'dataFile' path below if needed.
%   3. From MATLAB: run  nested_lasso_ms_fatigue
%
% The script prints performance metrics and writes:
%   - ElasticNet_Firth_Wald.csv : Firth refit coefficients, SEs, z, p, OR, CIs

clear; clc; rng(7,'twister');

%% 1) Load data -----------------------------------------------------------

% Relative or absolute path to the CSV file (edit for your machine)
dataFile = fullfile('data','lasso_with_caudateandSLF.csv');

if ~exist(dataFile, 'file')
    error('Data file not found: %s. Please check the path.', dataFile);
end

data = readtable(dataFile, 'VariableNamingRule','preserve');

% Outcome (assumes first column is the binary label)
outcomeVar = data.Properties.VariableNames{1};   % e.g., 'Dummy MS'
y_all      = data{:, outcomeVar};                % REQUIRED coding: 1=fatigued, 0=non-fatigued
assert(all(ismember(unique(y_all(~isnan(y_all))), [0 1])), ...
       'Outcome must be coded 0/1 (0=non-fatigued, 1=fatigued).');

% Candidate predictors (must match CSV headers exactly)
clin_vars = {'age','sex_num','disease_duration','education','MADRS', ...
             'MSQOL_54_PH','MSQOL_54_MH','9-HPT_right','NUM_LES','VOL_LES'};
vol_vars  = {'WM','CAU','THAL'};
tms_vars  = {'CSP','SICI','ICF','RMT','MEP latency'};
eeg_vars  = {'Bp2p_left_pow_avg','Bp2p_right_pow_avg','Bp2p_front_pow_avg'};
fa_vars   = {'FA_CST_R','FA_CST_L','FA_SLF_R','FA_SLF_L'};

predictorNames = [clin_vars, vol_vars, tms_vars, eeg_vars, fa_vars];

% Keep only predictors that exist in the file
existsMask = ismember(predictorNames, data.Properties.VariableNames);
if ~all(existsMask)
    warning('Missing predictors removed: %s', ...
            strjoin(predictorNames(~existsMask), ', '));
    predictorNames = predictorNames(existsMask);
end

% Design matrix
Xmat_all = data{:, predictorNames};

% ---- COMPLETE-CASE (listwise deletion), as in the primary analysis ------
rowKeep = all(isfinite(Xmat_all), 2) & isfinite(y_all);
Xmat    = Xmat_all(rowKeep, :);
y       = y_all(rowKeep);

fprintf('N included after listwise deletion: %d (of %d total)\n', numel(y), numel(y_all));
% -------------------------------------------------------------------------

%% 2) CV configuration ----------------------------------------------------

Kouter     = 5;                                    % outer folds
Kinner     = 5;                                    % inner folds for tuning
AlphaGrid  = [0.1 0.3 0.5 0.7 0.9 1];              % elastic-net to LASSO (α=1)
NumLambda  = 100;
useWeights = false;                                % set true for class weighting

%% 3) One full nested CV run (pooled held-out metrics + selection) --------

res = run_one_cv(Xmat, y, Kouter, Kinner, AlphaGrid, NumLambda, useWeights);

fprintf('\n=== Single nested-CV run (pooled outer folds) ===\n');
fprintf('Overall AUC: %.3f\n', res.AUC);
fprintf('PR-AUC: %.3f\n', res.PRAUC);
fprintf('Accuracy: %.3f, Balanced accuracy: %.3f\n', res.accuracy, res.bal_accuracy);
fprintf('Brier score: %.4f\n', res.Brier);

selFreq = res.selFreq(:);
tblSel  = table(predictorNames(:), selFreq, ...
                'VariableNames', {'Predictor','SelectionFreq'});
disp(tblSel);

%% 4) Repeated outer CV for CIs and stability -----------------------------

nRepeats = 50;                   % increase for tighter CIs
metrics  = zeros(nRepeats,5);    % [AUC PRAUC acc balacc Brier]
selMat   = zeros(nRepeats, numel(predictorNames));

rng(7,'twister');
for r = 1:nRepeats
    rr = run_one_cv(Xmat, y, Kouter, Kinner, AlphaGrid, NumLambda, useWeights);
    metrics(r,:) = [rr.AUC, rr.PRAUC, rr.accuracy, rr.bal_accuracy, rr.Brier];
    selMat(r,:)  = rr.selFreq(:)';   % selection freq for this run
end

avg = mean(metrics,1);
ci  = prctile(metrics,[2.5 97.5],1);   % 95% CI

fprintf('\n=== Repeated nested-CV (n=%d repeats) ===\n', nRepeats);
fprintf('AUC mean=%.3f [%.3f, %.3f]\n',  avg(1), ci(1,1), ci(2,1));
fprintf('PR-AUC mean=%.3f [%.3f, %.3f]\n',avg(2), ci(1,2), ci(2,2));
fprintf('Acc mean=%.3f  BalAcc mean=%.3f  Brier mean=%.3f\n', ...
        avg(3), avg(4), avg(5));

meanSel = mean(selMat,1)';  % average selection frequency across repeats
selTable = table(predictorNames(:), meanSel, ...
                 'VariableNames', {'Predictor','MeanSelFreq'});
selTable = sortrows(selTable, 'MeanSelFreq', 'descend');
disp(selTable);

%% 5) Descriptive penalized refit on ≥50% selection frequency ------------
% Outcome coding reminder:
%   y = 1 -> fatigued (MS_F)
%   y = 0 -> non-fatigued (MS_NF)
% So OR < 1 corresponds to lower odds of fatigue (potentially "protective").

% Choose selection vector:
% keep single-run frequencies if you wan
