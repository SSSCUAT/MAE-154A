%% ========================================================================
% DYNAMIC PLOTTER FOR MULTIPLE CSV ITERATION FILES
% ========================================================================
% 1) Find all CSV files in the current folder
files = dir('*.csv');

if isempty(files)
    error('No CSV files found in the current folder!');
end

% 2) Sort the files by their creation/modification date (Earliest first)
[~, sort_idx] = sort([files.datenum]);
files = files(sort_idx);

% Initialize empty arrays to hold data from ALL files combined
all_weights = [];
all_results = [];

% 3) Loop through each file chronologically, read it, and stack the data
for k = 1:length(files)
    filename = files(k).name;
    fprintf('Reading data from: %s\n', filename);
    
    % Read the current CSV
    current_data = readtable(filename);
    
    % Append this file's data to our master lists
    all_weights = [all_weights; current_data.W_total];
    
    % Safely append the results, regardless of how readtable formatted them
    % Convert to string array to standardize, then stack
    all_results = [all_results; string(current_data.Result)];
end

% 4) Create the iteration numbers based on the TOTAL number of runs across all files
total_runs = length(all_weights);
iterations = (1:total_runs)';

% 5) BULLETPROOF RESULT PARSING
% Check if the text is '1'/'0', 'true'/'false', or 'Pass'/'Fail'
pass_idx = (all_results == "1") | strcmpi(all_results, "true") | strcmpi(all_results, "Pass");
fail_idx = (all_results == "0") | strcmpi(all_results, "false") | strcmpi(all_results, "Fail");

% 6) Plotting Setup
figure('Name', 'Dynamic MonteCarlo');
hold on;
grid on;

% Plot Passes (Blue dots) and Fails (Red dots)
if any(pass_idx)
    scatter(iterations(pass_idx), all_weights(pass_idx), 2, 'b', 'filled');
end
if any(fail_idx)
    scatter(iterations(fail_idx), all_weights(fail_idx), 2, 'r', 'filled');
end

% Format the graph
xlabel('Total Iteration Number (Chronological)', 'FontWeight', 'bold');
ylabel('Total Aircraft Weight (lbs)', 'FontWeight', 'bold');
title(sprintf('Monte Carlo Results: %d Total Runs', total_runs));

% Only add to legend what actually exists
if any(pass_idx) && any(fail_idx)
    legend({'Pass', 'Fail'}, 'Location', 'best');
elseif any(pass_idx)
    legend({'Pass'}, 'Location', 'best');
else
    legend({'Fail'}, 'Location', 'best');
end

xlim([0, total_runs + 1]);
hold off;

% Print a quick summary to the command window
success_rate = (sum(pass_idx) / total_runs) * 100;
fprintf('-----------------------------------\n');
fprintf('Total Runs Processed: %d\n', total_runs);
fprintf('Total Passes: %d\n', sum(pass_idx));
fprintf('Success Rate: %.1f%%\n', success_rate);
fprintf('-----------------------------------\n');