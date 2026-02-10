function out = circleAimpoint(in)
% ============================================================
% Waypoint guidance and target intercept using Proportional Navigation (PN)
% Author: D. Toohey
%
% PURPOSE:
% This function computes:
%   - A 2D aimpoint (East, North)
%   - A commanded altitude
%   - A commanded velocity
%
% The UAV operates in three modes:
%   1) SEARCH    -> Follow predefined waypoints
%   2) INTERCEPT -> PN guidance toward target
%   3) CLIMB     -> Post-intercept climb maneuver
%
% MODE LOGIC SUMMARY:
% - If target is farther than DETECTION_RANGE:
%       -> SEARCH (waypoints)
% - If target enters DETECTION_RANGE:
%       -> INTERCEPT (PN guidance)
% - If target is captured (within CAPTURE_RADIUS):
%       -> CLIMB for fixed time
% - After climb:
%       -> Return to SEARCH
%
% NOTE:
% - Circling behavior has been removed
% - This is Simulink-safe and uses persistent memory
% ============================================================

%% ================= INPUTS =================
% Target position (inertial frame)
tar_E = in(1);        % Target East position
tar_N = in(2);        % Target North position

% UAV position
pE = in(3);           % UAV East position
pN = in(4);           % UAV North position

% Waypoint index (passed in externally)
way_num = in(5);      % Current waypoint number

% Simulation time
simTime = in(6);      % Simulation time [s]

% Target velocity
tar_VE = in(7);       % Target East velocity
tar_VN = in(8);       % Target North velocity

%% ================= PERSISTENT MEMORY =================
% These variables retain their values between function calls

persistent mode               % Guidance mode
persistent interceptCount     % Number of completed intercepts
persistent interceptedFlag    % Indicates recent intercept
persistent climbStartTime     % Time when climb mode started

% Initialize persistent variables on first call
if isempty(mode)
    mode = 1;                 % 1=SEARCH, 2=INTERCEPT, 3=CLIMB
    interceptCount = 0;       % No intercepts initially
    interceptedFlag = 0;      % Target not yet intercepted
    climbStartTime = 0;       % No climb started yet
end

%% ================= PARAMETERS =================
% Detection and capture thresholds %changed units to ft manually 
DETECTION_RANGE = 6000;       % Range to begin intercept [ft]
CAPTURE_RADIUS  = 200;         % Capture distance [ft]

% Velocity limits
V_CRUISE = 130;               % Cruise speed [ft/s]
V_MIN    = 80;                % Minimum allowed speed
V_MAX    = 180;               % Maximum allowed speed
V_SLOW_RADIUS = 50;          % Radius for slowing down near target

% Altitude commands
ALT_SEARCH    = 200;         % Altitude during waypoint search
ALT_INTERCEPT = 30;          % Altitude during intercept
ALT_CLIMB     = 100;         % Altitude during climb-out

% Climb behavior
CLIMB_TIME = 10;              % Time to remain in climb mode [s]

% Proportional Navigation parameters
N = 3;                        % Navigation constant
AIM_DIST = 1000;              % Distance to project aimpoint ahead

%% ================= GEOMETRY =================
% Relative position from UAV to target
rel_E = tar_E - pE;           % Relative East distance
rel_N = tar_N - pN;           % Relative North distance

% Slant range to target
range = sqrt(rel_E^2 + rel_N^2);

% Target ground speed magnitude
tarSpeed = sqrt(tar_VE^2 + tar_VN^2);

%% ================= WAYPOINTS =================
% Fixed waypoint list [East, North]
waypoints = [ ...
    1000   5000;
    5000  10000;
    10000 11000;
    15000 12000 ];

%% ================= MODE TRANSITIONS =================
% SEARCH -> INTERCEPT
if mode == 1 && range < DETECTION_RANGE &&  interceptCount == 0
   mode = 2;
end

% INTERCEPT -> CLIMB
if mode == 2 && range < CAPTURE_RADIUS
    mode = 3;
    interceptCount = interceptCount + 1;   % Log intercept
    interceptedFlag = 1;                   % Set intercept flag
    climbStartTime = simTime;               % Record climb start time
end

% CLIMB -> SEARCH (after fixed time)
if mode == 3 && simTime > climbStartTime + CLIMB_TIME
    mode = 1;               % back to SEARCH
    interceptedFlag = 0;    
    way_num = 1;            % reset to first waypoint to start patrolling again
end

%% ================= SEARCH MODE =================
if mode == 1
    % Current waypoint position
    tempE = waypoints(way_num,1);
    tempN = waypoints(way_num,2);

    % Distance from UAV to waypoint
    dist_wp = sqrt((pE - tempE)^2 + (pN - tempN)^2);

    % If waypoint reached, advance to next
    if dist_wp < 50
        way_num = way_num + 1;

        % Loop back to first waypoint
        if way_num > size(waypoints,1)
            way_num = 1;
            interceptCount = 0;%to fix the intercet mode not engaging on reruns
        end
    end

    % Output commands for SEARCH mode
    out(1) = tempE;           % Aimpoint East
    out(2) = tempN;           % Aimpoint North
    out(3) = way_num;         % Waypoint index
    out(4) = ALT_SEARCH;      % Commanded altitude
    out(5) = V_CRUISE;        % Commanded speed
    return
end

%% ================= INTERCEPT MODE =================
if mode == 2
    % Line-of-sight (LOS) angle to target
    lambda = atan2(rel_E, rel_N);

    % LOS rate (lambda_dot)
    if range < 10
        lambda_dot = 0;       % Avoid division by zero
    else
        lambda_dot = (rel_E * tar_VN - rel_N * tar_VE) / (range^2);
    end

    % PN heading command
    heading_cmd = lambda + N * lambda_dot;

    % Aimpoint projected forward along commanded heading
    aimpoint_E = pE + AIM_DIST * sin(heading_cmd);
    aimpoint_N = pN + AIM_DIST * cos(heading_cmd);
%% Available velocity control (wont implement until we finalize plane parameters)
    % Velocity control near target
    if range < V_SLOW_RADIUS
        vCommand = min(V_MAX, tarSpeed + 20); %max(V_MIN, V_CRUISE * range / V_SLOW_RADIUS);
    else
        vCommand = min(V_MAX, tarSpeed + 20);
    end
%%
    % Output commands for INTERCEPT mode
    out(1) = aimpoint_E;      % Aimpoint East
    out(2) = aimpoint_N;      % Aimpoint North
    out(3) = way_num;         % Waypoint index (unchanged)
    out(4) = ALT_INTERCEPT;   % Commanded altitude
    out(5) = vCommand;        % Commanded speed
    return
end

%% ================= CLIMB MODE =================
if mode == 3
    % Simple forward nudge during climb
    out(1) = pE + 100;        % Aimpoint East
    out(2) = pN + 100;        % Aimpoint North
    out(3) = way_num;         % Waypoint index
    out(4) = ALT_CLIMB;       % Commanded altitude
    out(5) = V_MAX;           % Max speed
    return
end

end


         
         
   
   
       
