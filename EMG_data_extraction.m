%% EMG ANALYSIS – MS TAPPING STUDY
%
% Estrazione e quantificazione dell’EMG (FDI) dai dati grezzi g.Tec (.hdf5)
%
% Background:
%   Durante la raccolta dati, in una parte dei soggetti l’ordine dei canali
%   era diverso a causa di un errore di collegamento (cable mistake) fra il
%   cappuccio EEG e l’amplificatore. Per questi soggetti (oldchanlocs) i
%   canali sono stati registrati con una mappa di elettrodi diversa.
%
%   In questa analisi EMG:
%     - I soggetti con cablaggio originale (“oldchanlocs”) usano il file:
%           chanlocs_68_gTEC
%     - I soggetti con cablaggio corretto (“new”) usano:
%           new_chanlocs_68_gTEC
%
%   In questo modo, le posizioni dei canali (inclusi EMG/ECG) sono
%   coerenti all’interno di ciascun gruppo e correttamente etichettate.
%
% Output per soggetto:
%   - AUC (per trial) del segnale EMG filtrato e rettificato
%   - Durata EMG (in secondi) per ogni trial
%   - Indici dei punti scelti (inizio/fine) per ogni trial
%
% Script originale: Elisa Tatti
% Ultima revisione: 2025

clear; close all;

%% Paths
eeglab_path = '/Users/etatti/City College Dropbox/Elisa Tatti/eeglab2024.2.1/';
addpath(genpath(eeglab_path));

folderext = '/Users/etatti/City College Dropbox/Elisa Tatti/Data_storage_Tapping_EEG_MS_Siena/TAPPING_EEG/Externalchans/';

% Scegli la/e condizioni da analizzare (1=MS_F, 2=MS_NF, 3=HC)
for c = 1   % cambia in [1 2 3] se vuoi tutte le condizioni

    switch c
        case 1
            cond       = '1';
            subjects   = [1:9 12:19];   % 18 soggetti, Sbj10 no EEG
            oldchanlocs = [1:9];        % soggetti con file chanlocs “vecchio”
            g          = 1;             % gruppo 1

        case 2
            cond       = '2';
            subjects   = 12;            % es. soggetto singolo per debug
            oldchanlocs = [1:14];
            g          = 2;

        case 3
            cond       = '3';
            subjects   = 11;            % es. soggetto singolo per debug
            oldchanlocs = [1 2 3:7];
            g          = 3;
    end

    %% Cartelle di output EMG
    folderemg  = ['/Users/etatti/City College Dropbox/Elisa Tatti/Data_storage_Tapping_EEG_MS_Siena/TAPPING_EEG/EMG_data/Group', cond, '/'];
    if ~exist(folderemg, 'dir');   mkdir(folderemg);   end

    foldertrial = [folderemg, 'Singletrial_EMG/'];
    if ~exist(foldertrial, 'dir'); mkdir(foldertrial); end

    %% Lista file nella directory corrente (deve essere quella con gli .hdf5)
    currentdir = pwd;
    files      = dir(currentdir);
    names      = {files.name}';

    % Ordina alfabeticamente
    [~, alorder]   = sort(lower(names));
    newstructure    = names(alorder);

    % Avvia EEGLAB
    eeglab;

    count = 1;   % contatore soggetti per AUC_avg_all (se lo usi)

    for s = subjects

        %% Cerca il file hdf5 corrispondente al soggetto s
        foundFile = '';
        for i = 1:numel(newstructure)
            if numel(newstructure{i}) > 20
                % qui ti basi sul fatto che il nome inizia con l’ID paziente (es “11_…”)
                if contains(newstructure{i}(1:6), num2str(s)) && ...
                        strcmp(newstructure{i}(end-4:end), '.hdf5')
                    foundFile = newstructure{i};
                    break;
                end
            end
        end

        if isempty(foundFile)
            fprintf('Nessun file hdf5 trovato per Sbj%d (cond %s)\n', s, cond);
            continue;
        end

        fprintf('\n=== Sbj%d, cond %s ===\n', s, cond);
        fprintf('File: %s\n', foundFile);

        %% Carica file .hdf5 in EEGLAB
        [ALLEEG, EEG, ~, ALLCOM] = eeglab; %#ok<ASGLU>
        EEG = pop_loadhdf5('filename', foundFile, ...
                           'filepath', currentdir, ...
                           'rejectchans', [], ...
                           'ref_ch',[]);
        eeglab redraw;

        %% Determina condizione “TAPPING_GroupX” in base a g
        substringtap = {'tapping','Tapping','TAPPING'};
        if any(contains({foundFile}, substringtap))
            tapCond = ['TAPPING_Group', num2str(g)];
            name    = ['Sbj', num2str(s), '_', tapCond];
        else
            error('Non riesco a determinare la condizione TAPPING per %s', foundFile);
        end

        %% Chanlocs: old vs new (cable mistake)
        % Se il soggetto è in oldchanlocs, usa chanlocs_68_gTEC (vecchio cablaggio)
        % altrimenti usa new_chanlocs_68_gTEC (cablaggio corretto).
        if ismember(s, oldchanlocs)
            expcond = 'old';
        else
            expcond = 'new';
        end

        if strcmp(expcond, 'new')
            EEG = pop_chanedit(EEG, ...
                'load', {'/Users/etatti/City College Dropbox/Elisa Tatti/Data_storage_Tapping_EEG_MS_Siena/chanlocs/new_chanlocs_68_gTEC', ...
                'filetype','xyz'}, ...
                'changefield', {64,'datachan',1}, ...
                'changefield', {64,'type',''});
        else
            EEG = pop_chanedit(EEG, ...
                'load', {'/Users/etatti/City College Dropbox/Elisa Tatti/Data_storage_Tapping_EEG_MS_Siena/chanlocs/chanlocs_68_gTEC', ...
                'filetype','xyz'}, ...
                'changefield', {29,'datachan',1}, ...
                'changefield', {29,'type',''});
        end

        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, 1);
        EEG = eeg_checkset(EEG);
        eeglab redraw;

        %% Eccezioni in cui Trigger 1 / 2 sono invertiti
        % Gruppo 1: Sbj11
        % Gruppo 2: Sbj12
        % Gruppo 3: Sbj7 e Sbj8
        is_swapped = (strcmp(tapCond,'TAPPING_Group1') && s == 11) || ...
                     (strcmp(tapCond,'TAPPING_Group2') && s == 12) || ...
                     (strcmp(tapCond,'TAPPING_Group3') && (s == 7 || s == 8));

        if is_swapped
            trgName = 'Trigger 2';
        else
            trgName = 'Trigger 1';
        end

        % Elimina i primi 10 trigger del tipo selezionato
        prova = find(strcmp({EEG.event.type}, trgName));
        if numel(prova) >= 10
            EEG.event(prova(1:10)) = [];
        end
        eeglab redraw;

        % Elimina l'ultimo trigger dello stesso tipo (fine recording)
        prova = find(strcmp({EEG.event.type}, trgName));
        if ~isempty(prova)
            EEG.event(prova(end)) = [];
        end
        eeglab redraw;

        %% Rimuovi i primi 2 s di recording
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0); %#ok<NASGU>
        EEG = eeg_checkset(EEG);
        EEG = pop_select(EEG, 'notime', [0 2]);  % rimuovi 0–2 s

        %% Seleziona canale EMG (EMG2 = FDI)
        % Canali esterni: 65–68 = ECG1 ECG2 EMG1 EMG2
        EEG = pop_select(EEG, 'channel', {'EMG2'});
        eeglab redraw;

        %% Filtro EMG
        % High-pass 20 Hz (locutoff), mantenendo le alte frequenze
        EEG = pop_eegfiltnew(EEG, 'locutoff', 20);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2, ...
            'setname', [name, '_filt'], 'gui', 'off');

        % Notch 49–51 Hz
        EEG = pop_eegfiltnew(EEG, 'locutoff', 49, 'hicutoff', 51, 'revfilt', 1);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 3, ...
            'setname', [name, '_filt_notch'], 'gui', 'off');
        EEG = eeg_checkset(EEG);

        %% Epoching rispetto al trigger corretto
        if is_swapped
            EEG = pop_epoch(EEG, {'Trigger 2'}, [0 1.5], ...
                'newname', [name, '_epochs'], 'epochinfo', 'yes');
        else
            EEG = pop_epoch(EEG, {'Trigger 1'}, [0 1.5], ...
                'newname', [name, '_epochs'], 'epochinfo', 'yes');
        end
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 4, ...
            'setname', [name, '_epochs'], 'gui', 'off');
        EEG = eeg_checkset(EEG);

        %% Estrai EMG [samples x trials]
        EMG = squeeze(EEG.data(:, :, :));   % [1 x time x trials] -> [time x trials]
        fs  = EEG.srate;
        nqf = fs / 2;
        time = (0:(size(EMG,1) - 1)) / fs;

        fco = 20;   % cutoff per l’envelope (low-pass 20 Hz dopo rettifica)

        AUC        = nan(size(EMG,2), 1);
        EMG_length = nan(size(EMG,2), 1);
        Trial_ind  = nan(size(EMG,2), 2);

        for t = 1:size(EMG, 2)
            fprintf('  Trial %d / %d\n', t, size(EMG,2));

            %% Rettifica EMG
            emg_rect = abs(EMG(:,t) - nanmean(EMG(:,t)));

            % Butterworth 2° ordine low-pass su rectified
            [b, a] = butter(2, fco * 1.25 / nqf);   % 1.25 = correzione per filtfilt
            z = filtfilt(b, a, emg_rect);           % linear envelope

            %% Soglia e ricerca automatica finestra
            bas  = z(1:size(EMG,1));     % qui hai usato tutto il trial come “baseline”
            meant = mean(bas);
            SDt   = 0.3 * meant;
            start = meant + SDt;

            % Finestra di inizio (es. 100–257)
            emgwindow_beg    = z(100:257);
            activity_indices = find(emgwindow_beg > start);

            if isempty(activity_indices)
                beginning_index = [];
            else
                beg_val         = emgwindow_beg(activity_indices(1));
                beginning_index = find(z == beg_val, 1, 'first');
            end

            % Se non trovi nulla, passi comunque al selezionamento manuale
            if isempty(beginning_index) || beginning_index + 27 >= numel(z)
                beginning_index = round(0.2 * fs);  % fallback a 200 ms
            end

            emgwindow_fin       = z(beginning_index+27 : min(beginning_index+27+309-100, numel(z)));
            activity_indices_end = find(emgwindow_fin < meant);

            if isempty(activity_indices_end)
                fin_index = beginning_index + round(0.2 * fs);   % fallback 200 ms dopo
            else
                fin_val   = emgwindow_fin(min(activity_indices_end));
                fin_index = find(z == fin_val, 1, 'first');
            end

            %% Plot automatica e conferma
            figure;
            plot(z); hold on;
            plot(beginning_index, z(beginning_index), 'ro', 'MarkerSize', 10);
            plot(fin_index,       z(fin_index),       'ko', 'MarkerSize', 10);
            title(sprintf('Sbj%d %s – Trial %d', s, tapCond, t));

            prompt = 'Are you okay with the selected points? (y/n): ';
            x = input(prompt, 's');

            if strcmpi(x, 'y')
                Trial_ind(t,:) = [beginning_index, fin_index];
                time_points    = beginning_index:fin_index;
                selected_signal = z(time_points);
                close(gcf);
            else
                % Selezione manuale con ginput
                while true
                    figure;
                    plot(z); hold on;
                    title('Click on the beginning and end points of the activity');
                    [xclick, ~] = ginput(2);

                    [~, begin_index] = min(abs((1:numel(z)) - xclick(1)));
                    [~, end_index]   = min(abs((1:numel(z)) - xclick(2)));

                    Trial_ind(t,:) = [begin_index, end_index];
                    time_points    = begin_index:end_index;
                    selected_signal = z(time_points);

                    hold on;
                    plot(begin_index, z(begin_index), 'go', 'MarkerSize', 10);
                    plot(end_index,   z(end_index),   'ko', 'MarkerSize', 10);

                    fprintf('Selected points (trial %d): [%d  %d]\n', t, begin_index, end_index);
                    answer = input('Are you okay with the selected points? (y/n): ', 's');

                    if strcmpi(answer, 'y')
                        close(gcf);
                        break;
                    else
                        close(gcf);
                        Trial_ind(t,:) = [NaN NaN];
                        time_points    = [];
                    end
                end
            end

            %% AUC e durata EMG per il trial
            if ~isempty(time_points)
                area_under_curve = trapz(time_points, selected_signal);
                AUC(t,:)         = area_under_curve;

                % durata EMG in secondi
                EMG_length(t,:)  = time(time_points(end)) - time(time_points(1));
            end

            clear area_under_curve time_points begin_index end_index emg_rect z ...
                  beginning_index fin_index emgwindow_fin activity_indices_end;
        end

        %% Salvataggio per soggetto (decommenta se vuoi salvare)
        % save([folderemg,'SingleTrials_Time_EMG_', name, '_EMG_FDI.mat'], 'EMG_length');
        % save([folderemg,'SingleTrials_AUC_',  name, '_EMG_FDI.mat'], 'AUC');
        % save([folderemg,'Idx_EMG_',          name, '_EMG_FDI.mat'], 'Trial_ind');

        % Media AUC per soggetto (se vuoi costruire un vettore gruppo)
        % sum_auc        = nansum(AUC, 1);
        % AUC_avg        = sum_auc / size(EMG, 2);
        % AUC_avg_all(count, :) = AUC_avg;
        % count = count + 1;

        % save([folderemg,'All_sbj_AUC_Group', num2str(g), '_EMG_FDI.mat'], 'AUC_avg_all');

        clear AUC EMG_length Trial_ind EMG;
    end
end
