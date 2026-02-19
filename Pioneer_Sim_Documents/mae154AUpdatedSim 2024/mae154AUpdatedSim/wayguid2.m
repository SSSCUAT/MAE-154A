function out = wayguid2(in)

tar_E = in(1);
tar_N = in(2);
pE = in(3);
pN = in(4);
heading = in(5);

delta_E = tar_E - pE;
delta_N = tar_N - pN;

tar_head = atan2(delta_E, delta_N);

delta_psi = tar_head - heading;

% Wrap angle to [-pi, pi]
if delta_psi > pi
    delta_psi = delta_psi - 2*pi;
elseif delta_psi < -pi
    delta_psi = delta_psi + 2*pi;
end

phi_comm = 0.6 * delta_psi;

% Bank angle saturation
if phi_comm > 30*pi/180
    phi_comm = 30*pi/180;
elseif phi_comm < -30*pi/180
    phi_comm = -30*pi/180;
end

out(1) = phi_comm;
out(2) = tar_E;
out(3) = tar_N;

end

       
