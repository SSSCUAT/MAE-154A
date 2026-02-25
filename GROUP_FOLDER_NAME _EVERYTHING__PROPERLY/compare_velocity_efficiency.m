clear; clc; close all;

%% --- 1. USER SETTINGS ---
% APC files contain data for many RPMs. Set this to match one in your file.
% Popular values are 3000, 4000, 5000, 6000.
targetRPM = 5000; 

% This line tells MATLAB to look for all text files in the current folder.
fileList = dir('*.txt'); 

% Error check: make sure you actually have files in the folder!
if isempty(fileList)
    error('No .txt files found! Ensure your propeller files are in this folder.');
end

%% --- 2. PLOT PREPARATION ---
% Create a high-quality figure window
figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
hold on; grid on; box on;
colors = lines(length(fileList)); % Assigns a different color to each file

fprintf('--- ANALYSIS RESULTS FOR %d RPM ---\n', targetRPM);

%% --- 3. THE EXTRACTION LOOP ---
% We will loop through every file found in the folder one by one.
for i = 1:length(fileList)
    fileName = fileList(i).name;
    fid = fopen(fileName, 'r'); % 'r' means open the file for "Reading"
    
    data = [];           % This will store our final table of numbers
    isCorrectRPM = false; % This is a "switch" to tell us when to start saving data
    
    % The 'while' loop reads the file line-by-line until it reaches the end (feof)
    while ~feof(fid)
        currentLine = fgetl(fid); 
        
        % STEP A: Look for the RPM header in the text file
        if contains(currentLine, 'PROP RPM =')
            val = sscanf(currentLine, ' PROP RPM = %f');
            if abs(val - targetRPM) < 10 % If the file RPM matches our target
                isCorrectRPM = true;      % Turn our "switch" ON
                fgetl(fid); fgetl(fid);  % Skip the 2 header lines with column names
                continue;                % Go to the next line to start getting numbers
            else
                isCorrectRPM = false;     % Turn "switch" OFF for other RPM sections
            end
        end
        
        % STEP B: If the switch is ON, we grab the columns of numbers
        if isCorrectRPM
            % APC data layout: [V(mph), J, Ct, Cp, Efficiency, ...]
            numRow = sscanf(currentLine, '%f %f %f %f %f');
            
            % If the line is empty or contains text, the RPM data block has ended
            if isempty(numRow)
                isCorrectRPM = false; 
            else
                data = [data; numRow']; % Add the new row of numbers to our matrix
            end
        end
    end
    fclose(fid); % Close the file after reading it
    
    %% --- 4. CALCULATIONS & PLOTTING ---
    if ~isempty(data)
        % Extract specific columns from our data matrix
        v_mph = data(:, 1);      % Velocity in Miles per Hour
        v_ms  = v_mph * 0.44704; % Convert mph to Meters per Second (Engineering Standard)
        efficiency = data(:, 5);  % Propulsive Efficiency (0.0 to 1.0)
        
        % Find the "Sweet Spot": Max Efficiency and the Speed it occurs at
        [maxEff, idx] = max(efficiency);
        bestV = v_ms(idx);
        
        % Output the "Analysis" results to the MATLAB command window
        fprintf('File: %s\n', fileName);
        fprintf('  > Max Efficiency: %.2f%%\n', maxEff * 100);
        fprintf('  > Optimal Cruise Speed: %.2f m/s\n\n', bestV);
        
        % Create the efficiency curve plot
        plot(v_ms, efficiency, 'LineWidth', 2.5, 'Color', colors(i,:), ...
             'DisplayName', strrep(fileName, '.txt', ''));
        
        % Place a circle on the peak efficiency point
        plot(bestV, maxEff, 'ok', 'MarkerFaceColor', colors(i,:), 'HandleVisibility', 'off');
    end
end

%% --- 5. GRAPH FORMATTING ---
% These lines make the graph look professional for your project report
xlabel('Flight Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Efficiency (\eta)', 'FontSize', 12, 'FontWeight', 'bold');
title(['UAV Propeller Efficiency Comparison at ', num2str(targetRPM), ' RPM'], 'FontSize', 14);

% Add a legend using the filenames (underscores removed for neatness)
legend('Location', 'northeastoutside', 'Interpreter', 'none');

% Set Y-axis from 0 to 1 (Efficiency can't be more than 100%)
ylim([0 1]); 
grid minor;