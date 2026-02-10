%waypoint guidance
%D.Toohey

function  out = wayguid2(in)

tar_E = in(1);
tar_N = in(2);
pE = in(3);
pN = in(4);
heading = in(5);

% if way_num == 0
%     way_num = 1;
% end

% %pre-defined waypoints  (pE   pN
% waypoints = [0 5000;...
%              -1000 10000;...
%              5000 10000;...
%              10000 10000;...
%              12000 8000;...
%              10000 6000;...
%              8000  3000];
%          
% distance threshold used to switch to next waypoint
% dist_thresh = 50;
         
% tar_E = waypoints(way_num,1);
% tar_N = waypoints(way_num,2);

delta_E = tar_E - pE;
delta_N = tar_N - pN;
% way_dist = (delta_E^2 + delta_N^2)^.5;
% if way_dist < dist_thresh
%     way_num = way_num + 1;
% end

tar_head = atan2(delta_E,delta_N);

delta_psi = tar_head - heading;


%check for angles larger than 180 deg
if delta_psi > pi
    delta_psi = delta_psi - 2*pi;
elseif delta_psi < -pi
    delta_psi = delta_psi + 2*pi;
end

phi_comm = .6*delta_psi;

if phi_comm > 30*pi/180
    phi_comm = 30*pi/180;
elseif phi_comm < -30*pi/180;
    phi_comm = -30*pi/180;
end

%out(1) = way_num;
out(1) = phi_comm;
out(2) = tar_E;
out(3) = tar_N;
         
         
   
   
       
