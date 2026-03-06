%% ========================================================================
% SCRIPT 3: VISUAL SENSITIVITY ANALYSIS & PASSING OVERVIEW
% (Run this from INSIDE your Analysis folder)
% ========================================================================
clear; clc; close all;

% 1) Pop up a UI window to select WHICH extracted CSV to analyze
[file, path] = uigetfile('*.csv', 'Select the Extracted Passing Cases CSV to Analyze');

if isequal(file, 0)
    disp('Analysis canceled. No file selected.');
    return;
end

% Load the data
filename = fullfile(path, file);
fprintf('Loading data from: %s\n', file);
passing_data = readtable(filename);

% Extract just the name of the file for our plot titles
[~, name_only, ~] = fileparts(file); 

% =========================================================================
% GENERATE THE OVERVIEW PLOT OF ALL PASSING CASES
% =========================================================================
figure('Name', sprintf('%s - Overview', name_only));
hold on; grid on;

% Check if we saved the Original Iteration number, otherwise use a generic index
if ismember('Original_Iteration', passing_data.Properties.VariableNames)
    x_vals = passing_data.Original_Iteration;
    x_label = 'Original Iteration Number';
else
    x_vals = 1:height(passing_data);
    x_label = 'Passing Case Index';
end

% Plot the passes
scatter(x_vals, passing_data.W_total, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.7);

% Format the graph
xlabel(x_label, 'FontWeight', 'bold');
ylabel('Total Aircraft Weight (lbs)', 'FontWeight', 'bold');
title(sprintf('%s: W_total Overview (%d Passing Runs)', name_only, height(passing_data)), 'Interpreter', 'none');

% Auto-scale the Y-axis to perfectly frame the passing weights
y_margin = 0.5;
if max(passing_data.W_total) > min(passing_data.W_total)
    ylim([min(passing_data.W_total) - y_margin, max(passing_data.W_total) + y_margin]);
end
hold off;

% =========================================================================
% GENERATE THE INDIVIDUAL VARIABLE SCATTER PLOTS
% =========================================================================
vars_to_plot = {'L_fuse', 'x_wle', 'Arw', 'Art', 'bw'}; 

% -> THIS IS LIKELY WHERE YOUR ERROR STARTED (Line 58)
for i = 1:length(vars_to_plot)
    var_name = vars_to_plot{i};
    
    if ~ismember(var_name, passing_data.Properties.VariableNames)
        fprintf('Warning: Column "%s" not found. Skipping...\n', var_name);
        continue;
    end
    
    x_data = passing_data.(var_name);
    y_data = passing_data.W_total;
    
    figure('Name', sprintf('Weight vs %s', var_name));
    scatter(x_data, y_data, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    hold on; grid on;
    
    % Draw trendline
    if length(x_data) > 1 && max(x_data) > min(x_data)
        coeffs = polyfit(x_data, y_data, 1);
        x_fit = linspace(min(x_data), max(x_data), 100);
        y_fit = polyval(coeffs, x_fit);
        plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
        legend('Passing Designs', 'Trendline', 'Location', 'best');
    end
    
    xlabel(var_name, 'Interpreter', 'none', 'FontWeight', 'bold');
    ylabel('Total Weight (lbs)', 'FontWeight', 'bold');
    title(sprintf('%s: W_total vs. %s', name_only, var_name), 'Interpreter', 'none');
    
    x_margin = 0.1 * (max(x_data) - min(x_data));
    if x_margin == 0
        x_margin = 1; 
    end 
    xlim([min(x_data) - x_margin, max(x_data) + x_margin]);
    hold off;
end % -> THIS IS THE 'END' YOU WERE MISSING!

fprintf('Analysis complete for %s! Review plots to adjust your bounds.\n', file);