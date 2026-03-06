%% ========================================================================
% SCRIPT 1: RAW DATA PLOTTER & PASSING OVERVIEW (CUSTOM MERGE ORDER)
% ========================================================================
clear; clc; close all;

% 1) Let the user select files ONE BY ONE to guarantee the exact plot order
fileList = {};
pathList = {};
count = 1;

disp('Please select your files in the exact order you want them plotted.');
disp('Click "Cancel" when you are finished selecting files.');

while true
    prompt_str = sprintf('Select File #%d (or click Cancel to finish and plot)', count);
    [file, path] = uigetfile('*.csv', prompt_str);
    
    if isequal(file, 0)
        break; % User clicked cancel to finish selection
    end
    
    fileList{count} = file;
    pathList{count} = path;
    count = count + 1;
end

if isempty(fileList)
    disp('Plotting canceled. No files selected.');
    return;
end

% Print the order so the user can verify
fprintf('\nMerging %d files in this EXACT order:\n', length(fileList));
for i = 1:length(fileList)
    fprintf('  Plot Chunk %d: %s\n', i, fileList{i});
end

% 2) Load and stack the data in the user's chosen order
all_data = table();
for k = 1:length(fileList)
    filename = fullfile(pathList{k}, fileList{k});
    current_data = readtable(filename);
    
    if isempty(all_data)
        all_data = current_data;
    else
        common_cols = intersect(all_data.Properties.VariableNames, current_data.Properties.VariableNames, 'stable');
        all_data = [all_data(:, common_cols); current_data(:, common_cols)];
    end
end

total_runs = height(all_data);
all_data.Original_Iteration = (1:total_runs)';

% Parse passes and fails
results_str = string(all_data.Result);
pass_idx = (results_str == "1") | strcmpi(results_str, "true") | strcmpi(results_str, "Pass");
fail_idx = (results_str == "0") | strcmpi(results_str, "false") | strcmpi(results_str, "Fail");

% =========================================================================
% FIGURE 1: RAW DATA (PASSES AND FAILS)
% =========================================================================
figure('Name', 'Raw Monte Carlo Plot (All Runs)');
hold on; grid on;
if any(pass_idx)
    scatter(all_data.Original_Iteration(pass_idx), all_data.W_total(pass_idx), 2, 'b', 'filled');
end
if any(fail_idx)
    scatter(all_data.Original_Iteration(fail_idx), all_data.W_total(fail_idx), 2, 'r', 'filled');
end
xlabel('Total Iteration Number', 'FontWeight', 'bold');
ylabel('Total Aircraft Weight (lbs)', 'FontWeight', 'bold');
title(sprintf('Monte Carlo Results: %d Total Runs', total_runs));
legend({'Pass', 'Fail'}, 'Location', 'best');
xlim([0, total_runs + 1]);
hold off;

% =========================================================================
% FIGURE 2: ZOOMED-IN PASSES ONLY (ALL MERGED FILES)
% =========================================================================
if any(pass_idx)
    figure('Name', 'All Passing Cases Overview');
    hold on; grid on;
    
    % Extract just the passing data
    passing_data = all_data(pass_idx, :);
    
    % Plot passes with slightly larger markers
    scatter(passing_data.Original_Iteration, passing_data.W_total, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
    
    xlabel('Total Iteration Number', 'FontWeight', 'bold');
    ylabel('Total Aircraft Weight (lbs)', 'FontWeight', 'bold');
    title(sprintf('All Passing Cases: %d Total Passes', height(passing_data)));
    
    % Auto-scale the Y-axis to tightly hug the passing weights
    y_margin = 0.5;
    if max(passing_data.W_total) > min(passing_data.W_total)
        ylim([min(passing_data.W_total) - y_margin, max(passing_data.W_total) + y_margin]);
    end
    xlim([0, total_runs + 1]);
    hold off;
end

% =========================================================================
% CONSOLE SUMMARY & TOP 5 DATA FIGURE (FOR REPORT)
% =========================================================================
fprintf('\nOVERALL SUMMARY\n');
fprintf('-----------------------------------\n');
fprintf('Total Runs Processed: %d\n', total_runs);
fprintf('Total Passes: %d\n', sum(pass_idx));
fprintf('Success Rate: %.1f%%\n', (sum(pass_idx) / total_runs) * 100);

if any(pass_idx)
    % Sort passing data by total weight, lowest to highest
    passing_data_sorted = sortrows(passing_data, 'W_total');
    
    % Determine if we have at least 5 to show, or less
    num_to_display = min(5, height(passing_data_sorted));
    
    % Extract the top cases
    top_5_cases = passing_data_sorted(1:num_to_display, :);
    
    % --- CREATE THE POP-UP TABLE FIGURE ---
    fig_table = figure('Name', 'Top 5 Lowest Weight Passing Cases', 'Position', [100, 100, 1000, 200]);
    
    % Generate the UI Table inside the figure
    uitable(fig_table, 'Data', table2cell(top_5_cases), ...
        'ColumnName', top_5_cases.Properties.VariableNames, ...
        'RowName', 'numbered', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.05 0.96 0.90]); 
        
    fprintf('\nSUCCESS: Generated a pop-up table figure with the Top %d cases for your report!\n', num_to_display);
else
    fprintf('\nNo passing cases available to display.\n');
end
fprintf('-----------------------------------\n');