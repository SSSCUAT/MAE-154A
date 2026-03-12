function out = circleAimpoint(in)
% ============================================================
% Waypoint guidance and target intercept using Proportional Navigation (PN)
% Author: D. Toohey / Upgraded by You
% ============================================================

    % 1. Pre-allocate the output vector so Simulink C++ coder doesn't crash
    out = zeros(1, 5);

    %% ================= INPUTS =================
    tar_E = in(1);        % Target East position
    tar_N = in(2);        % Target North position
    pE = in(3);           % UAV East position
    pN = in(4);           % UAV North position
    way_num = in(5);      % Current waypoint number
    simTime = in(6);      % Simulation time [s]
    tar_VE = in(7);       % Target East velocity
    tar_VN = in(8);       % Target North velocity

    %% ================= ORIGINAL SAFETY OVERRIDE =================
    % Preserved from the original code. If the C++ visualization expects 
    % the sim to "park" at 500 seconds, this ensures it doesn't break.
    if simTime > 500
        out(1) = 0;
        out(2) = 0;
        out(3) = way_num;
        out(4) = 500;
        out(5) = 130;
        return; % Safe to return here because 'out' is already pre-allocated
    end

    %% ================= PERSISTENT MEMORY =================
    persistent mode               
    persistent interceptCount     
    persistent interceptedFlag    
    persistent climbStartTime     
    
    if isempty(mode)
        mode = 1;                 % 1=SEARCH, 2=INTERCEPT, 3=CLIMB
        interceptCount = 0;       
        interceptedFlag = 0;      
        climbStartTime = 0;       
    end

    %% ================= PARAMETERS =================
    DETECTION_RANGE = 6000;       % Range to begin intercept [ft]
    CAPTURE_RADIUS  = 200;        % Capture distance [ft]
    
    V_CRUISE = 130;               % Cruise speed [ft/s]
    V_MIN    = 80;                % Minimum allowed speed
    V_MAX    = 180;               % Maximum allowed speed
    V_SLOW_RADIUS = 50;           % Radius for slowing down near target
    
    ALT_SEARCH    = 200;          % Altitude during waypoint search
    ALT_INTERCEPT = 30;           % Altitude during intercept
    ALT_CLIMB     = 100;          % Altitude during climb-out
    
    CLIMB_TIME = 10;              % Time to remain in climb mode [s]
    N = 3;                        % Navigation constant
    AIM_DIST = 1000;              % Distance to project aimpoint ahead

    %% ================= GEOMETRY =================
    rel_E = tar_E - pE;           
    rel_N = tar_N - pN;           
    range = sqrt(rel_E^2 + rel_N^2);
    tarSpeed = sqrt(tar_VE^2 + tar_VN^2);

    %% ================= WAYPOINTS =================
    waypoints = [ ...
        1000   5000;
        5000  10000;
        10000 11000;
        15000 12000 ];

    %% ================= MODE TRANSITIONS =================
    if mode == 1 && range < DETECTION_RANGE && interceptCount == 0
       mode = 2;
    end
    
    if mode == 2 && range < CAPTURE_RADIUS
        mode = 3;
        interceptCount = interceptCount + 1;   
        interceptedFlag = 1;                   
        climbStartTime = simTime;               
    end
    
    if mode == 3 && simTime > climbStartTime + CLIMB_TIME
        mode = 1;               
        interceptedFlag = 0;    
        way_num = 1;            
    end

    %% ================= GUIDANCE LOGIC (If/Elseif structure) =================
    % Using if/elseif guarantees Simulink assigns the output correctly 
    % without getting confused by early 'return' statements.

    if mode == 1
        % --- SEARCH MODE ---
        tempE = waypoints(way_num,1);
        tempN = waypoints(way_num,2);
        
        dist_wp = sqrt((pE - tempE)^2 + (pN - tempN)^2);
        
        if dist_wp < 50
            way_num = way_num + 1;
            if way_num > size(waypoints,1)
                way_num = 1;
                interceptCount = 0; 
            end
        end
        
        out(1) = tempE;           
        out(2) = tempN;           
        out(3) = way_num;         
        out(4) = ALT_SEARCH;      
        out(5) = V_CRUISE;        

    elseif mode == 2
        % --- INTERCEPT MODE (Predicted Lead Fix) ---
        
        % Approximate time to reach the target (Time-to-go)
        % Adding a small number (0.1) prevents division by zero if V_CRUISE drops
        t_go = range / (V_CRUISE + 0.1); 
        
        % Lead Pursuit: Project where the target will be based on its velocity.
        % We use 'N' as an aggressiveness gain. If N=3, it aggressively leads the target.
        % N=1 is a pure true-intercept prediction.
        lead_factor = N / 3; 
        
        aimpoint_E = tar_E + (tar_VE * t_go * lead_factor);
        aimpoint_N = tar_N + (tar_VN * t_go * lead_factor);

        % Velocity control near target
        if range < V_SLOW_RADIUS
            vCommand = min(V_MAX, tarSpeed + 20); 
        else
            vCommand = min(V_MAX, tarSpeed + 20);
        end

        out(1) = aimpoint_E;      
        out(2) = aimpoint_N;      
        out(3) = way_num;         
        out(4) = ALT_INTERCEPT;   
        out(5) = vCommand;      

    elseif mode == 3
        % --- CLIMB MODE ---
        out(1) = pE + 100;        
        out(2) = pN + 100;        
        out(3) = way_num;         
        out(4) = ALT_CLIMB;       
        out(5) = V_MAX;           
    end
end
