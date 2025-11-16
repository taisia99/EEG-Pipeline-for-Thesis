%% STEP: Clean Events and Save Continuous Data with Target Triggers
clear; clc;

%% === Paths ===
input_paths  = {...
    '/Users/Taisia/Documents/lab_work/MPI/EEG/ENCprocessed_data/after_ica', ...
    '/Users/Taisia/Documents/lab_work/MPI/EEG/ENCprocessed_data/interpolated_data'};
output_path = '/Users/Taisia/Documents/lab_work/MPI/EEG/ENCprocessed_data/target_triggers'; % Change this path if needed

if ~exist(output_path, 'dir')
    mkdir(output_path);
end

%% === Parameters ===
file_ext_pattern = '*_afterICA.set';

%% === Collect files ===
files = [];
for p = 1:length(input_paths)
    files_tmp = dir(fullfile(input_paths{p}, file_ext_pattern));
    if ~isempty(files_tmp)
        files = files_tmp;
        input_path = input_paths{p};
        break;
    end
end
if isempty(files)
    error('No files found in either after_ica or interpolated_data folders!');
end

%% === Main loop across participants ===
for subject = 1:length(files)
    EEG = pop_loadset('filename', files(subject).name, 'filepath', input_path);
    [~, subject_id, ~] = fileparts(files(subject).name);
    fprintf('\n--- Processing %s ---\n', subject_id);

    %% === Event cleaning ===
    if isfield(EEG, 'event') && ~isempty(EEG.event)

        % Convert all event types to strings for consistent processing
        all_types = string({EEG.event.type});

        % First, remove unwanted or invalid event markers
        bad_labels = ["999", "boundary", "actiCAP Data On"];
        good_idx = ~ismember(all_types, bad_labels);
        EEG.event = EEG.event(good_idx);
        EEG = eeg_checkset(EEG);

        % Keep only stimulus events of the form 0_* or 1_*
        all_types = string({EEG.event.type}); % update after removal
        stim_idx = startsWith(all_types, "0_") | startsWith(all_types, "1_");
        EEG.event = EEG.event(stim_idx);
        EEG = eeg_checkset(EEG);

        fprintf('Events after filtering: %d\n', length(EEG.event));

        % Display remaining events (useful for debugging)
        disp(EEG.event)

        % Latency correction (optional)
        % Uncomment or modify if latency transformation is required
        %
        % latencies_corrected = round([EEG.event.latency] / 100000 * 250);
        % for i = 1:length(EEG.event)
        %     EEG.event(i).latency = latencies_corrected(i);
        % end

    else
        warning('EEG.event is empty or invalid for %s.', subject_id);
    end

    %% === Save continuous data with filtered target triggers ===
    if isfield(EEG, 'event') && ~isempty(EEG.event)
        % Generate output filename
        save_name = strrep(files(subject).name, '_afterICA.set', '_target_triggers.set');
        EEG = pop_saveset(EEG, 'filename', save_name, 'filepath', output_path);

        fprintf('Saved continuous data with target triggers: %s\n', ...
            fullfile(output_path, save_name));
    else
        warning('Nothing to save — EEG.event is empty after filtering!');
    end

end
