% map display function -- D.Toohey
function out = show_map(in)
    pn = in(1);
    pe = in(2);
    alt = -in(3);
    tar_E = in(4);
    tar_N = in(5);
    padre_E = in(6); % Previously whale_E
    padre_N = in(7); % Previously whale_N
    
    % Persistent variables to maintain data between time steps
    persistent initialized
    % Figure 1 vars (Clean View)
    persistent hUAV1 hTar1 hPadre1 hText intercepts is_intercepting stored_waypoints
    % Figure 2 vars (Trail View)
    persistent hUAV2 hPadre2 hTar2 uav_hist_E uav_hist_N padre_hist_E padre_hist_N
    
    % FIX: We now check if the handle is empty OR if it's no longer a valid graphic object
    if isempty(initialized) || isempty(hUAV1) || ~isvalid(hUAV1)
        initialized = true;
        brown = [0.6 0.3 0]; % Padres brown
        
        % ==========================================
        % FIGURE 1: Clean View (No Trails)
        % ==========================================
        figure(1);
        clf(1);
        hold on;
        
        stored_waypoints = [tar_E, tar_N];
        
        hUAV1 = plot(pe, pn, '.b', 'MarkerSize', 16);
        hTar1 = plot(stored_waypoints(:,1), stored_waypoints(:,2), '.r', 'MarkerSize', 16);
        hPadre1 = plot(padre_E, padre_N, '.', 'Color', brown, 'MarkerSize', 20);
        
        title('Figure 1: Live Tracking (No Trails)');
        xlabel('pE');
        ylabel('pN');
        ylim([-15000 15000]);
        xlim([-15000 15000]);
        axis equal;
        grid on;
        legend({'UAV position', 'Target Waypoints', 'Padres Fan'}, 'Location', 'northwest', 'AutoUpdate', 'off');
        
        intercepts = 0;
        is_intercepting = false; 
        hText = text(-14000, 13000, 'Interceptions: 0', 'FontSize', 12, 'FontWeight', 'bold');
        hold off;
        
        % ==========================================
        % FIGURE 2: Full Path History (Trails)
        % ==========================================
        figure(2);
        clf(2);
        hold on;
        
        % Initialize history arrays
        uav_hist_E = pe;
        uav_hist_N = pn;
        padre_hist_E = padre_E;
        padre_hist_N = padre_N;
        
        % Plot the first points
        hUAV2 = plot(uav_hist_E, uav_hist_N, '.b');
        hTar2 = plot(tar_E, tar_N, '.r', 'MarkerSize', 16);
        hPadre2 = plot(padre_hist_E, padre_hist_N, '.', 'Color', brown);
        
        title('Figure 2: Full Path History');
        xlabel('pE');
        ylabel('pN');
        ylim([-15000 15000]);
        xlim([-15000 15000]);
        axis equal;
        grid on;
        legend({'UAV path', 'Target position', 'Padres Fan path'}, 'Location', 'northwest', 'AutoUpdate', 'off');
        hold off;
        
    else
        % ==========================================
        % UPDATE FIGURE 1 (Clean View)
        % ==========================================
        set(hUAV1, 'XData', pe, 'YData', pn);
        set(hPadre1, 'XData', padre_E, 'YData', padre_N);
        
        % Update waypoints (keep up to 4)
        last_wp = stored_waypoints(end, :);
        if (tar_E ~= last_wp(1) || tar_N ~= last_wp(2)) && size(stored_waypoints, 1) < 4
            stored_waypoints = [stored_waypoints; tar_E, tar_N];
            set(hTar1, 'XData', stored_waypoints(:,1), 'YData', stored_waypoints(:,2));
        end
        
        % Interception logic (500ft capture radius)
        dist = sqrt((pe - padre_E)^2 + (pn - padre_N)^2);
        if dist < 500
            if ~is_intercepting
                intercepts = intercepts + 1;
                is_intercepting = true;
                set(hText, 'String', sprintf('Interceptions: %d', intercepts));
            end
        else
            is_intercepting = false; 
        end
        
        % ==========================================
        % UPDATE FIGURE 2 (Full History)
        % ==========================================
        % Add new positions to the history arrays
        uav_hist_E = [uav_hist_E, pe];
        uav_hist_N = [uav_hist_N, pn];
        padre_hist_E = [padre_hist_E, padre_E];
        padre_hist_N = [padre_hist_N, padre_N];
        
        % Update the plots with the full arrays to create the trail
        set(hUAV2, 'XData', uav_hist_E, 'YData', uav_hist_N);
        set(hPadre2, 'XData', padre_hist_E, 'YData', padre_hist_N);
        set(hTar2, 'XData', tar_E, 'YData', tar_N); % Just shows current target
    end
    
    % Slower speed limit to match original feel
    pause(0.05); 
    drawnow;
    
    out = pn;
end