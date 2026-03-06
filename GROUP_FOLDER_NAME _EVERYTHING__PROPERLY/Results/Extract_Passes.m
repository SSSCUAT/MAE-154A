%% ========================================================================
% SCRIPT 2: EXTRACT & SAVE PASSING CASES
% ========================================================================
clear; clc;

% 1) Pop up a UI to let the user select specific raw CSV files
[fileList, path] = uigetfile('*.csv', 'Select the Raw Monte Carlo CSV files', 'MultiSelect', 'on');

if isequal(fileList, 0)
    disp('Extraction canceled. No files selected.');
    return;
end

% Convert single file selection to cell array to prevent loop errors
if ischar(fileList)
    fileList = {fileList};
end

% 2) Load and stack ONLY the selected files
all_data = table();
fprintf('Extracting data from %d selected files...\n', length(fileList));

for k = 1:length(fileList)
    filename = fullfile(path, fileList{k});
    current_data = readtable(filename);
    
    if isempty(all_data)
        all_data = current_data;
    else
        common_cols = intersect(all_data.Properties.VariableNames, current_data.Properties.VariableNames, 'stable');
        all_data = [all_data(:, common_cols); current_data(:, common_cols)];
    end
end

% 3) Filter for passes
results_str = string(all_data.Result);
pass_idx = (results_str == "1") | strcmpi(results_str, "true") | strcmpi(results_str, "Pass");

if ~any(pass_idx)
    fprintf('No passing cases found in the selected files!\n');
    return;
end

passing_data = all_data(pass_idx, :);

% 4) Ask the user what to name this extraction (e.g., "Iter_1")
prompt = {'Enter a unique name for this run (e.g., Iter_1):'};
dlgtitle = 'Name Your Extracted Data';
run_name_input = inputdlg(prompt, dlgtitle, [1 40], {'Iter_1'});

if isempty(run_name_input)
    disp('Save canceled by user.');
    return;
end

run_name = run_name_input{1};

% 5) Save the passing cases into the Analysis folder
output_folder = fullfile(pwd, 'Analysis');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

output_filename = fullfile(output_folder, sprintf('Passing_Cases_%s.csv', run_name));
writetable(passing_data, output_filename);

fprintf('SUCCESS: Extracted %d passing cases and saved to:\n%s\n', height(passing_data), output_filename);