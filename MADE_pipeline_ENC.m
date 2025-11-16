%% EEG Preprocessing Script
clear; clc;

%% === Initialize Parameters ===
adjust_time_offset = 0; % Initialize if not set (0 = NO, 1 = YES)
filter_timeoffset = 0; % anti-aliasing filter time offset (ms)
stimulus_timeoffset = 0; % stimulus time offset (ms)
response_timeoffset = 0; % response time offset (ms)
stimulus_markers = {}; % stimulus markers
response_markers = {}; % response markers
down_sample = 1; % 0 = NO, 1 = YES
sampling_rate = 250; % target sampling rate
save_interim_result = 1; % 0 = NO, 1 = YES
output_format = 1; % 1 = .set, 2 = .mat
rerefer_data = 1; % 0 = NO, 1 = YES
reref = []; % reference channels (empty for average reference)
voltthres_rejection = 1; % 0 = NO, 1 = YES
volt_threshold = [-100 100]; % [lower upper] threshold (μV)
interp_epoch = 1; % 0 = NO, 1 = YES
interp_channels = 1; % 0 = NO, 1 = YES
remove_baseline = 1; % 0 = NO, 1 = YES
baseline_window = [-100 0]; % baseline window (ms)
output_location = '/Users/Taisia/Documents/lab_work/MPI/EEG/processed_data';

%% === Path Setup ===
eeglab_path = '/Users/Taisia/Documents/lab_work/MPI/EEG/eeglab2020_0';
addpath(genpath(eeglab_path));
eeglab; close;

rawdata_location = '/Users/Taisia/Documents/lab_work/MPI/EEG';
data_root = fullfile(rawdata_location, 'mydata');

main_chanlocs_path = fullfile(eeglab_path, 'sample_locs', 'MEME_32_REF.bvef');
extra_chanlocs_path = fullfile(eeglab_path, 'sample_locs', 'MEME_32_REF_A1A2.bvef');

if ~isfile(main_chanlocs_path); error(['Main channel locations file not found: ' main_chanlocs_path]); end
if ~isfile(extra_chanlocs_path); error(['Extra channel locations file not found: ' extra_chanlocs_path]); end

chanlocs_main = loadbvef(main_chanlocs_path);
chanlocs_extra = loadbvef(extra_chanlocs_path);

%% === Load and Process EEG Data Files ===
file_list = dir(fullfile(data_root, '**', 'sub_*_ENC.vhdr'));

disp(['Number of files found: ', num2str(length(file_list))]);

datafile_names = fullfile({file_list.folder}, {file_list.name});

%% Initialize output variables
reference_used_for_faster = cell(length(datafile_names), 1);
faster_bad_channels = cell(length(datafile_names), 1);
ica_preparation_bad_channels = cell(length(datafile_names), 1);
length_ica_data = zeros(length(datafile_names), 1);
total_ICs = zeros(length(datafile_names), 1);
ICs_removed = cell(length(datafile_names), 1);
total_epochs_before_artifact_rejection = zeros(length(datafile_names), 1);
total_epochs_after_artifact_rejection = zeros(length(datafile_names), 1);
total_channels_interpolated = zeros(length(datafile_names), 1);

%% Create output folders
if save_interim_result == 1
    if ~exist(fullfile(output_location, 'filtered_data'), 'dir')
        mkdir(fullfile(output_location, 'filtered_data'));
    end
    if ~exist(fullfile(output_location, 'ica_data'), 'dir')
        mkdir(fullfile(output_location, 'ica_data'));
    end
end
if ~exist(fullfile(output_location, 'processed_data'), 'dir')
    mkdir(fullfile(output_location, 'processed_data'));
end

%% Loop over all data files
for subject = 1:length(datafile_names)
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names{subject});
    
    %% STEP 1: Import EEG data (RAW)
    try
        EEG = pop_loadbv([], datafile_names{subject});
        EEG = eeg_checkset(EEG);
        
        % === Save raw immediately after import ===
        if save_interim_result == 1
            [~, name, ~] = fileparts(datafile_names{subject});
            pop_saveset(EEG, 'filename', [name '_raw.set'], ...
                        'filepath', fullfile(output_location, 'filtered_data'));
        end
    catch
        error('Failed to load file: %s', datafile_names{subject});
    end
    
    %% STEP 2: Import and merge channel locations
    % Display channel information
    disp('--- Main channel location file (chanlocs_main):');
    disp({chanlocs_main.labels});
    
    disp('--- Extra channel location file (chanlocs_extra):');
    disp({chanlocs_extra.labels});
    
    % Merge channels without duplicates
    all_labels = {chanlocs_main.labels};
    for k = 1:length(chanlocs_extra)
        if ~ismember(chanlocs_extra(k).labels, all_labels)
            chanlocs_main(end+1) = chanlocs_extra(k);
            all_labels{end+1} = chanlocs_extra(k).labels;
        end
    end
    
    % Create placeholder labels if EEG.chanlocs is empty
    if isempty({EEG.chanlocs.labels}) || all(cellfun(@isempty, {EEG.chanlocs.labels}))
        EEG.chanlocs = struct('labels', arrayfun(@(i) ['chan' num2str(i)], 1:size(EEG.data,1), 'UniformOutput', false));
        warning('EEG.chanlocs.labels were empty - created placeholders "chan1", "chan2", ...');
    end
    
    % Filter only channels present in EEG
    data_labels = {EEG.chanlocs.labels};
    chanlocs_filtered = chanlocs_main([]); % empty structure with correct format
    
    for j = 1:length(chanlocs_main)
        if ismember(chanlocs_main(j).labels, data_labels)
            chanlocs_filtered(end+1) = chanlocs_main(j);
        end
    end
    
    % Check channel count match
    if length(chanlocs_filtered) ~= size(EEG.data,1)
        error('Error: Number of channels in EEG.data (%d) does not match filtered locations (%d)', ...
            size(EEG.data,1), length(chanlocs_filtered));
    else
        fprintf('Success: Assigned %d channels with coordinates.\n', length(chanlocs_filtered));
    end
    
    % Assign and check
    EEG.chanlocs = chanlocs_filtered;
    EEG = eeg_checkset(EEG);
    
    % Save after channel assignment
    if save_interim_result == 1
        [~, name, ~] = fileparts(datafile_names{subject});
        pop_saveset(EEG, 'filename', [name '_chanlocs.set'], 'filepath', fullfile(output_location, 'filtered_data'));
    end
    
    %% STEP 2.5: Relabel events, clean unwanted types, and assign stim1-stim96

% 1. Ensure all EEG.event(ev).type are strings without extra spaces
for ev = 1:length(EEG.event)
    if isnumeric(EEG.event(ev).type)
        EEG.event(ev).type = num2str(EEG.event(ev).type);
    end
    EEG.event(ev).type = strtrim(EEG.event(ev).type);
end

% 2. Display unique event types before cleaning
uniqueTypesBefore = unique({EEG.event.type});
fprintf('Found %d unique event types (before cleaning):\n', numel(uniqueTypesBefore));
disp(uniqueTypesBefore);

% 3. Remove unwanted events: 'S  3', 'S  4', 'S  5', 'boundary'
unwanted = {'S  3', 'S  4', 'S  5', 'boundary'};
keep_idx = ~ismember({EEG.event.type}, unwanted);
EEG.event = EEG.event(keep_idx);

% 4. Replace all 'S 12' with stim1...stim96
idx_S12 = find(strcmp({EEG.event.type}, 'S 12'));
fprintf('Found %d events with type ''S 12''.\n', numel(idx_S12));
if length(idx_S12) ~= 96
    error('Expected 96 ''S 12'' events, found %d.', length(idx_S12));
end
for k = 1:96
    EEG.event(idx_S12(k)).type = sprintf('stim%d', k);
end

% 5. Cheking types
uniqueTypesAfter = unique({EEG.event.type});
fprintf('Events after relabeling:\n');
disp(uniqueTypesAfter);

    
    %% STEP 4: Change sampling rate
    if down_sample == 1
        if floor(sampling_rate) > EEG.srate
            error('Sampling rate cannot be higher than recorded sampling rate');
        elseif floor(sampling_rate) ~= EEG.srate
            EEG = pop_resample(EEG, sampling_rate);
            EEG = eeg_checkset(EEG);
            
            % Save after downsampling
            if save_interim_result == 1
                [~, name, ~] = fileparts(datafile_names{subject});
                pop_saveset(EEG, 'filename', [name '_downsampled.set'], 'filepath', fullfile(output_location, 'filtered_data'));
            end
        end
    end
    
 %% STEP 6: Filter data
% Filtering
highpass = 0.3;
lowpass = 30;
high_transband = 1;
low_transband = 10;


hp_fl_order = ceil((3.3 / (high_transband / EEG.srate)) / 2) * 2;
lp_fl_order = ceil((3.3 / (low_transband / EEG.srate)) / 2) * 2;
hp_fl_order = max(hp_fl_order, 20);
lp_fl_order = max(lp_fl_order, 350);


high_cutoff = highpass + high_transband/2;
low_cutoff = lowpass + low_transband/2;

%  high-pass
EEG = eeg_checkset(EEG);
EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass', ...
                'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
EEG = eeg_checkset(EEG);

%  low-pass
EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass', ...
                'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);

% Notch filter 50 Гц
EEG = pop_eegfiltnew(EEG, 'locutoff', 49, 'hicutoff', 51, ...
                     'revfilt', 1, 'plotfreqz', 0);
EEG = eeg_checkset(EEG);

% Saving
if save_interim_result == 1
    [~, name, ~] = fileparts(datafile_names{subject});
    pop_saveset(EEG, 'filename', [name '_filtered.set'], ...
                      'filepath', fullfile(output_location, 'filtered_data'));
end
disp(['Current events at STEP 6: ' strjoin(unique({EEG.event.type}), ', ')]);

%% STEP 7: Rereference
if rerefer_data == 1
    EEG = eeg_checkset(EEG);
    EEG = pop_reref(EEG, reref);  %  = average reference
    EEG.etc.reref_info = struct(...
        'method', 'average', ...
        'original_ref', reref, ...
        'date', datestr(now));
    
    if save_interim_result == 1
        [~, name, ~] = fileparts(datafile_names{subject});
        pop_saveset(EEG, 'filename', [name '_reref.set'], ...
                          'filepath', fullfile(output_location, 'filtered_data'));
    end
end

%% Remove predefined bad channels before ICA

% for W7DZFN_T1 channels_to_remove = { 'FC6','Fp1', 'A2', 'P3'};
% for LF4RXR channels_to_remove = {'Fp1', 'T7','T8', 'P4', 'Pz'};
% for A546TQ channels_to_remove = {'Fp1', 'T7','CP5','CP6' };
% for CNPTG5: channels_to_remove = {'Fp1', 'Cz','T7','T8', 'F4', 'O1', 'Fp2',  'A2', 'A1', 'FT9'};
% for 3FZK2X: channels_to_remove = {'Fp1', 'FC6','T7', 'CP2' };
%for 3FZK2X: 
channels_to_remove = {'Fp1', 'FT10','T8', 'CP2' };

EEG_full = EEG;  % Save original EEG with all channels

remove_idx = find(ismember({EEG.chanlocs.labels}, channels_to_remove));

if isempty(remove_idx)
    warning('No channels marked for removal.');
    removed_channel_labels = {};
else
    if any(remove_idx > length(EEG.chanlocs)) || any(remove_idx < 1)
        error('Invalid channel index found.');
    end

    % Store removed channel labels for later interpolation
    removed_channel_labels = {EEG.chanlocs(remove_idx).labels};

    fprintf('Removed %d predefined channels: %s\n', ...
        numel(remove_idx), strjoin(removed_channel_labels, ', '));

    EEG = pop_select(EEG, 'nochannel', removed_channel_labels);
    EEG = eeg_checkset(EEG);
end

% %% STEP 7: Enhanced channel quality detection

EEG_copy = eeg_checkset(EEG);


    % Define frontal channels
    frontal_channels = {'Fz','F3','F4','F7','F8'};
    frontal_idx = find(ismember({EEG.chanlocs.labels}, frontal_channels));
    other_idx = setdiff(1:EEG.nbchan, frontal_idx);
    
    % Threshold values
    thresholds = struct(...
        'frontal', struct('variance', 4.0, 'corr', -2.5, 'hf_noise', 4.0), ...
        'other',   struct('variance', 3.0, 'corr', -3.0, 'hf_noise', 3.0) ...
    );
    
    % Calculate quality metrics
    [spectra, freqs] = spectopo(EEG_copy.data, 0, EEG_copy.srate, 'plot', 'off');
    metrics.variance = var(EEG_copy.data, [], 2);
    
    try
        C = corrcoef(EEG_copy.data');
        C(isnan(C)) = 0;
        metrics.corr = median(C, 2);
    catch
        warning('Could not calculate channel correlation. Replacing with zeros.');
        metrics.corr = zeros(EEG.nbchan, 1);
    end
    
    metrics.hf_noise = mean(spectra(:, freqs > 20 & freqs < 40), 2);
    
    % Find bad channels
    bad_chans = false(EEG.nbchan, 1);
    bad_metrics_count = zeros(EEG.nbchan, 1);
    
    for m = fieldnames(thresholds.other)'
        z = zscore(metrics.(m{1}));
        bad_chans(other_idx) = bad_chans(other_idx) | (abs(z(other_idx)) > thresholds.other.(m{1}));
        bad_metrics_count(other_idx) = bad_metrics_count(other_idx) + (abs(z(other_idx)) > thresholds.other.(m{1}));
    end
    
    for m = fieldnames(thresholds.frontal)'
        z = zscore(metrics.(m{1}));
        bad_chans(frontal_idx) = bad_chans(frontal_idx) | (abs(z(frontal_idx)) > thresholds.frontal.(m{1}));
        bad_metrics_count(frontal_idx) = bad_metrics_count(frontal_idx) + (abs(z(frontal_idx)) > thresholds.frontal.(m{1}));
    end
    
    % Additional correlation check (non-frontal)
    corr_bad = find(metrics.corr < median(metrics.corr) - 3*std(metrics.corr));
    bad_chans(setdiff(corr_bad, frontal_idx)) = true;
    
    % Decision to remove (non-frontal and ≥2 bad metrics)
    to_remove = (bad_metrics_count >= 2);
    to_remove(frontal_idx) = false;
    
    % Too many bad channels - enable "gentle" mode
    if sum(bad_chans) > EEG.nbchan * 0.8
        warning('Too many bad channels - skipping removal, interpolating noisy channels only');
        to_remove = false(size(to_remove)); % Don't remove channels
    end
    
    % Separate lists
    noisy_channels = find(bad_chans);
    removed_channels = find(to_remove);
    interpolated_channels = setdiff(noisy_channels, removed_channels);
    
    noisy_labels = {EEG.chanlocs(noisy_channels).labels};
    removed_labels = {EEG.chanlocs(removed_channels).labels};
    interp_labels = {EEG.chanlocs(interpolated_channels).labels};
    
    % Visualization
    if ~isempty(noisy_labels)
        eegplot(EEG.data, 'srate', EEG.srate, 'winlength', 10, ...
            'title', sprintf('Noisy channels\nRemoved: %s\nInterpolated: %s', ...
            strjoin(removed_labels, ', '), strjoin(interp_labels, ', ')));
    end
    
    % Only remove if not too many channels
too_many_bad = numel(removed_labels) > EEG.nbchan * 0.8;

if ~isempty(removed_labels) && ~too_many_bad
    fprintf('🗑️ Removed %d non-frontal channels: %s\n', ...
        length(removed_labels), strjoin(removed_labels, ', '));
    EEG = pop_select(EEG, 'nochannel', removed_labels);
    EEG = eeg_checkset(EEG);
elseif too_many_bad
    warning('Too many channels marked bad (%d of %d) - skipping removal, interpolating only.', ...
        numel(removed_labels), EEG.nbchan);
    removed_labels = {};  % Reset
end

% Check channel coordinates
has_coords = isfield(EEG.chanlocs, 'X') && isfield(EEG.chanlocs, 'Y') && isfield(EEG.chanlocs, 'Z');
has_nan_coords = false;
if has_coords
    coords = [[EEG.chanlocs.X]; [EEG.chanlocs.Y]; [EEG.chanlocs.Z]];
    has_nan_coords = any(isnan(coords(:))) || any(abs(coords(:)) > 1e3);
end

    
  
    %% STEP 9: Prepare data for ICA
    EEG_copy = eeg_checkset(EEG);
    
    % Perform 1Hz high pass filter
    transband = 1;
    fl_cutoff = highpass + (high_transband / 2);
    fl_order = 3.3 / (transband / EEG.srate);
    
    if mod(floor(fl_order),2) == 0
        fl_order = floor(fl_order);
    else
        fl_order = floor(fl_order) + 1;
    end
    
    EEG_copy = pop_firws(EEG_copy, 'fcutoff', fl_cutoff, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', fl_order, 'minphase', 0);
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Create 1 second epochs
    EEG_copy = eeg_regepochs(EEG_copy, 'recurrence', 1, 'limits', [0 1], 'rmbase', [NaN], 'eventtype', '999');
    EEG_copy = eeg_checkset(EEG_copy);
    
    % Find bad epochs
    vol_thrs = [-1000 1000]; % [lower upper] threshold limit(s) in uV.
    emg_thrs = [-100 30]; % [lower upper] threshold limit(s) in dB.
    emg_freqs_limit = [20 40]; % [lower upper] frequency limit(s) in Hz.
    
    % Find channels with >20% artifacted epochs
    chanCounter = 1; 
    ica_prep_badChans = [];
    numEpochs = EEG_copy.trials;
    all_bad_channels = 0;
    
    for ch = 1:EEG_copy.nbchan
        % Detect outlier voltage
        EEG_copy = pop_eegthresh(EEG_copy, 1, ch, vol_thrs(1), vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax, 0, 0);
        EEG_copy = eeg_checkset(EEG_copy);
        
        % Detect EMG artifacts
        EEG_copy = pop_rejspec(EEG_copy, 1, 'elecrange', ch, 'method', 'fft', 'threshold', emg_thrs, 'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);
        
        disp(class(EEG));

        % Find number of artifacted epochs
        EEG_copy = eeg_checkset(EEG_copy);
        EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
        artifacted_epochs = EEG_copy.reject.rejglobal;
        
        % Mark bad channels (>20% bad epochs)
        if sum(artifacted_epochs) > (numEpochs*20/100)
            ica_prep_badChans(chanCounter) = ch;
            chanCounter = chanCounter + 1;
        end
    end
    
    % Handle case where all channels are bad
    if numel(ica_prep_badChans) == EEG.nbchan || numel(ica_prep_badChans)+1 == EEG.nbchan
        all_bad_channels = 1;
        warning(['No usable data for datafile ', datafile_names{subject}]);
        
        % Save dataset with bad channels info
        if output_format == 1
            EEG = eeg_checkset(EEG);
            [~, name, ~] = fileparts(datafile_names{subject});
            EEG = pop_editset(EEG, 'setname', [name, '_no_usable_data_all_bad_channels']);
            EEG = pop_saveset(EEG, 'filename', [name, '_no_usable_data_all_bad_channels.set'], 'filepath', fullfile(output_location, 'processed_data'));
        else
            [~, name, ~] = fileparts(datafile_names{subject});
            save(fullfile(output_location, 'processed_data', [name, '_no_usable_data_all_bad_channels.mat']), 'EEG');
        end
       disp(['Current events at STEP 9: ' strjoin(unique({EEG.event.type}), ', ')])

        % Update output variables
        ica_preparation_bad_channels{subject} = num2str(1:EEG.nbchan);
        length_ica_data(subject) = 0;
        total_ICs(subject) = 0;
        ICs_removed{subject} = '0';
        total_epochs_before_artifact_rejection(subject) = 0;
        total_epochs_after_artifact_rejection(subject) = 0;
        total_channels_interpolated(subject) = 0;
        continue; % Skip to next subject
    else
        % Reject bad channels
        if ~isempty(ica_prep_badChans)
            EEG_copy = pop_select(EEG_copy, 'nochannel', ica_prep_badChans);
            EEG_copy = eeg_checkset(EEG_copy);
        end
    end
    
    % Update output variable
    if isempty(ica_prep_badChans)
        ica_preparation_bad_channels{subject} = '0';
    else
        ica_preparation_bad_channels{subject} = num2str(ica_prep_badChans);
    end
    
    % Reject artifacted epochs before ICA
    EEG_copy = pop_eegthresh(EEG_copy, 1, 1:EEG_copy.nbchan, vol_thrs(1), vol_thrs(2), EEG_copy.xmin, EEG_copy.xmax, 0, 0);
    EEG_copy = eeg_checkset(EEG_copy);
    
    EEG_copy = pop_rejspec(EEG_copy, 1, 'elecrange', 1:EEG_copy.nbchan, 'method', 'fft', 'threshold', emg_thrs, 'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);
    
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
    reject_artifacted_epochs = EEG_copy.reject.rejglobal;
    EEG_copy = pop_rejepoch(EEG_copy, reject_artifacted_epochs, 0);
    
    
 [~, name, ~] = fileparts(datafile_names{subject});
ica_file_path = fullfile(output_location, 'ica_data', [name, '_ica_data.set']);

% === ica already exists? ===
if output_format == 1 && exist(ica_file_path, 'file') == 2
    fprintf('ICA already exists for %s. Loading...\n', name);
    EEG_copy = pop_loadset('filename', [name, '_ica_data.set'], 'filepath', fullfile(output_location, 'ica_data'));
elseif output_format ~= 1 && exist(ica_mat_path, 'file') == 2
    fprintf('ICA already exists for %s. Loading...\n', name);
    load(ica_mat_path, 'EEG_copy');
else
    %% STEP 10: Run ICA
    length_ica_data(subject) = EEG_copy.trials;
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = pop_runica(EEG_copy, 'icatype', 'runica', 'extended', 1, 'stop', 1E-7, 'interupt', 'off');

    % Save ICA data if requested
    if save_interim_result == 1
        if output_format == 1
            EEG_copy = pop_editset(EEG_copy, 'setname', [name, '_ica_data']);
            EEG_copy = pop_saveset(EEG_copy, 'filename', [name, '_ica_data.set'], 'filepath', fullfile(output_location, 'ica_data'));
        else
            save(fullfile(output_location, 'ica_data', [name, '_ica_data.mat']), 'EEG_copy');
        end
    end
end

% ===transfer ICA components===
ICA_WINV = EEG_copy.icawinv;
ICA_SPHERE = EEG_copy.icasphere;
ICA_WEIGHTS = EEG_copy.icaweights;
ICA_CHANSIND = EEG_copy.icachansind;

EEG.icawinv = ICA_WINV;
EEG.icasphere = ICA_SPHERE;
EEG.icaweights = ICA_WEIGHTS;
EEG.icachansind = ICA_CHANSIND;
EEG = eeg_checkset(EEG);


total_ICs(subject) = size(EEG.icawinv, 2);

    
    %% STEP 11: Run ADJUST to detect artifact components
    try
        % Use modified ADJUST with relaxed thresholds
        [art, horiz, vert, blink, disc,...
            soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
            soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=ADJUST_mod(EEG, 'report', 'off');
        
        % Only mark components that are clearly artifacts
        ICs_removed{subject} = unique([blink; horiz; vert; disc]);
        
        % If too many components are marked (>30%), be more conservative
        if length(ICs_removed{subject}) > 0.3*total_ICs(subject)
            ICs_removed{subject} = unique([blink; horiz]); % Only remove eye artifacts
            warning('Too many components marked - only removing clear eye artifacts');
        end
    catch
        warning('ADJUST failed - using alternate artifact detection');
        % Alternative method using ICLabel
        EEG = iclabel(EEG);
        brain_thresh = 0.7;
        other_thresh = 0.8;
        
        brain_idx = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) < brain_thresh);
        other_idx = [];
        for i = 2:size(EEG.etc.ic_classification.ICLabel.classifications,2)
            other_idx = [other_idx; find(EEG.etc.ic_classification.ICLabel.classifications(:,i) > other_thresh)];
        end
        ICs_removed{subject} = unique(other_idx);
    end
    
   % Remove artifact components if any found
if ~isempty(ICs_removed{subject})
    EEG = pop_subcomp(EEG, ICs_removed{subject}, 0);
    EEG = eeg_checkset(EEG);

    % Save after ICA component removal
    if save_interim_result == 1
        [~, name, ~] = fileparts(datafile_names{subject});
        pop_saveset(EEG, 'filename', [name '_icaclean.set'], 'filepath', fullfile(output_location, 'filtered_data'));
    end
else
    ICs_removed{subject} = '0';
end


%% STEP 11: Interpolate previously removed channels after ICA
if exist('removed_channel_labels', 'var') && ~isempty(removed_channel_labels)
    % Load the saved full EEG with all original channels
    original_EEG = EEG_full;

    % Find indices of removed channels in original EEG
    interp_idx = find(ismember({original_EEG.chanlocs.labels}, removed_channel_labels));

    if isempty(interp_idx)
        warning('No matching channels found to interpolate.');
    else
        % Prepare EEG for interpolation
        EEG_interp = EEG;

        % Add empty channels (NaNs)
        n_interp = length(interp_idx);
        missing_data = nan(n_interp, size(EEG_interp.data, 2));
        EEG_interp.data = [EEG_interp.data; missing_data];
        EEG_interp.nbchan = size(EEG_interp.data, 1);

        % Add back chanlocs
        for i = 1:n_interp
            EEG_interp.chanlocs(end+1) = original_EEG.chanlocs(interp_idx(i));
        end

        % Reorder channels to match original
        [~, order] = ismember({original_EEG.chanlocs.labels}, {EEG_interp.chanlocs.labels});
        EEG_interp.data = EEG_interp.data(order, :);
        EEG_interp.chanlocs = EEG_interp.chanlocs(order);
        EEG_interp = eeg_checkset(EEG_interp);

        % Interpolate using spherical method
        EEG = eeg_interp(EEG_interp, interp_idx, 'spherical');
        EEG = eeg_checkset(EEG);

        fprintf('Interpolated %d channels: %s\n', ...
            length(interp_idx), strjoin({original_EEG.chanlocs(interp_idx).labels}, ', '));
    end
end

if save_interim_result == 1
    interpolated_folder = fullfile(output_location, 'interpolated_data');
    if ~exist(interpolated_folder, 'dir')
        mkdir(interpolated_folder);
    end

    [~, name, ~] = fileparts(datafile_names{subject});
    output_filename = [name '_after_ica_interp.set'];

    pop_saveset(EEG, 'filename', output_filename, 'filepath', interpolated_folder);
    fprintf('Saved dataset with interpolated channels: %s/%s\n', interpolated_folder, output_filename);
end


%% Preprocessing


% 1. : ±150 мкВ 
EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -150, 150, EEG.xmin, EEG.xmax, 0, 0); 

% 2. 
EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 5, 5, 0, 0); 

% 3. 
EEG = pop_rejkurt(EEG, 1, 1:EEG.nbchan, 5, 5, 0, 0);

% 4. 
EEG = pop_rejspec(EEG, 1, ...
    'elecrange', 1:EEG.nbchan, ...
    'method', 'fft', ...
    'threshold', [-50 50], ...      % could be [-35 35]
    'freqlimits', [1 50], ...
    'eegplotplotallrej', 0, ...
    'eegplotreject', 0);

% 5. 
EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);

% 6. 
bad_samples = find(EEG.reject.rejglobal);
if ~isempty(bad_samples)
    fprintf('Found %d bad segments\n', length(bad_samples));
    
end


% 
subject_id = 'sub_1VAEL7_ENC';

fprintf('\n--- Processing stim epochs ---\n');

% 'stim' with stim1, stim2 ...
event_types = {EEG.event.type};
pattern = '^stim\d+$';
matched_events = unique(event_types(~cellfun(@isempty, regexp(event_types, pattern))));

if isempty(matched_events)
    warning('No stim events found.');
else
    % creating epochs
    epoch_window = [-0.2, 1.0];

    try
        EEG_ep = pop_epoch(EEG, matched_events, epoch_window, 'newname', [subject_id '_stim'], 'epochinfo', 'yes');
        EEG_ep = eeg_checkset(EEG_ep);
        fprintf('Created %d stim epochs.\n', EEG_ep.trials);

        % Path
        raw_output_dir = fullfile(output_location, 'raw_epoched_data');
        if ~exist(raw_output_dir, 'dir')
            mkdir(raw_output_dir);
        end

        % saving file
        raw_filename = sprintf('%s_stim_raw.set', subject_id);
        pop_saveset(EEG_ep, 'filename', raw_filename, 'filepath', raw_output_dir);
        fprintf('Saved raw stim epochs to %s\n', fullfile(raw_output_dir, raw_filename));

    catch ME
        warning('Error creating stim epochs: %s', ME.message);
    end
end


end
