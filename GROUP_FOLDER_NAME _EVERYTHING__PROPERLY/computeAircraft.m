function outputs = computeAircraft(inputs)

%% =============================
% UNPACK INPUTS
% =============================
rho = inputs.rho;
mu  = inputs.mu;
a   = inputs.a;

W = inputs.W;

Arw     = inputs.Arw;
bw      = inputs.bw;
Sw      = inputs.Sw;
lambdaw = inputs.lambdaw;

Sth      = inputs.Sth;
bth      = inputs.bth;
lambdath = inputs.lambdath;
Art      = inputs.Art;

e        = inputs.e;
cla      = inputs.cla;
downwash = inputs.downwash;
it       = inputs.it;
Cmacw    = inputs.Cmacw;
tau      = inputs.tau;
CL_MAX   = inputs.CL_MAX;

x_wle    = inputs.x_wle;
x_cg_str = inputs.x_cg_str;
x_cg_tot = inputs.x_cg_tot;
x_cg_0   = inputs.x_cg_0;

Power = inputs.Power;
EF    = inputs.EF;
L_fuse = inputs.L_fuse;
W_fuse = inputs.W_fuse;
H_fuse = inputs.H_fuse;
%% (10) DISPLAY OPTIONS (directly from inputs)
makePlots            = inputs.makePlots;           
makeTable            = inputs.makeTable;           
makePrint_stability  = inputs.makePrint_stability; 
makePrint_tail_Volume = inputs.makePrint_tail_Volume;
%% Fuselage
Swet_f = 2*(L_fuse*W_fuse + L_fuse*H_fuse + W_fuse*H_fuse);
d_fuse = W_fuse;  % for fuselage form factor

%% ========================================================================
% 2) PLANFORM GEOMETRY (ROOT, TIP, MAC)
Crw = (2*Sw)/(bw*(1+lambdaw)); 
Ctw = lambdaw*Crw;
Cw  = (2/3)*Crw*((1+lambdaw+lambdaw^2)/(1+lambdaw));

Crth = (2*Sth)/(bth*(1+lambdath)); 
Ctth = lambdath*Crth;
Cth  = (2/3)*Crth*((1+lambdath+lambdath^2)/(1+lambdath));

%% ========================================================================
% 3) 3D LIFT CURVE SLOPES
aw = cla / (1 + (cla/(pi*e*Arw)));
at = cla / (1 + (cla/(pi*e*Art)));

%% ========================================================================
% 4) LONGITUDINAL REFERENCE LOCATIONS
hacw = 0.25; 
x_ac_w_LE = hacw*Cw;
Lh = 4*Cw; 
x_ac_t_LE = x_ac_w_LE + Lh; 
h_ac_t = x_ac_t_LE / Cw;

%% ========================================================================
% 5) NEUTRAL POINT + STATIC MARGIN
hn = (hacw + h_ac_t*((Sth/Sw)*(at/aw))*(1-downwash)) / (1 + (Sth/Sw)*(at/aw)*(1-downwash));
hcg_str = (x_cg_str - x_wle)/Cw; 
hcg_tot = (x_cg_tot - x_wle)/Cw; 
hcg_0 = (x_cg_0 - x_wle)/Cw;
SM_str = hn - hcg_str; 
SM_tot = hn - hcg_tot; 
SM_0 = hn - hcg_0;

%% ========================================================================
% 6) TAIL ARM + TAIL VOLUME
x_wac = x_wle + x_ac_w_LE; 
x_tac = x_wac + Lh;
lt_str = x_tac - x_cg_str; 
lt_tot = x_tac - x_cg_tot; 
lt_0 = x_tac - x_cg_0;
VH_str = (Sth*lt_str)/(Sw*Cw); 
VH_tot = (Sth*lt_tot)/(Sw*Cw); 
VH_0 = (Sth*lt_0)/(Sw*Cw);

%% ========================================================================
% 7) LONGITUDINAL COEFFICIENTS
VH  = VH_tot; 
hcg = hcg_tot;
CLalpha = aw + at*(Sth/Sw)*(1-downwash);
CL0     = -at*(Sth/Sw)*it;
Cmalpha = -CLalpha*(hn-hcg); 
Cm0 = Cmacw + VH*at*it;
CLdeltae = at*(Sth/Sw)*tau; 
Cmdeltae = -at*VH*tau;

%% ========================================================================
% 8) VELOCITY SWEEP
V = linspace(30, 300, 500); 
q = 0.5*rho*V.^2; 
M = V./a;

CL = W ./ (q.*Sw); 
alpha = CL ./ aw; 
alpha_t = (1-downwash).*alpha + it;
CLw = aw.*alpha; 
CLt = at.*alpha_t;

Re_w = rho.*V.*Cw ./ mu; 
Re_t = rho.*V.*Cth ./ mu; 
Re_f = rho.*V.*L_fuse ./ mu;

Swet_w = 2*Sw; 
Swet_t = 2*Sth;
Cf_w = 0.455 ./ ((log10(Re_w)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_t = 0.455 ./ ((log10(Re_t)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_f = 0.455 ./ ((log10(Re_f)).^2.58 .* (1 + 0.144*M.^2).^0.65);

tc = 0.12; xcm = 0.3; Lam = 0; 
K_surf = (1+(0.6/xcm)*tc+100*tc^4).*(1.34*M.^0.18 .* cos(Lam).^0.28);
Q_w = 1.0; Q_t = 1.05; Q_f = 1.0; 
K_f = 0.9 + 5/d_fuse^1.5 + d_fuse/400;

CD_misc = 0.002; 
CD_LP = 0.002;

CDp = K_surf.*Q_w.*Cf_w.*(Swet_w/Sw) + K_surf.*Q_t.*Cf_t.*(Swet_t/Sw) + K_f.*Q_f.*Cf_f.*(Swet_f/Sw) + CD_misc + CD_LP;
CDi_w = (CLw.^2)/(pi*Arw*e); 
CDi_t = (Sth/Sw).*(CLt.^2)/(pi*Art*e);
CDi_total = CDi_w + CDi_t; 
CD_total = CDp + CDi_total;

Dp = q.*Sw.*CDp; 
Di = q.*Sw.*CDi_total; 
D = q.*Sw.*CD_total;
L = q.*Sw.*CL; 
LD = L./D;

Power_Required_hp = (D.*V)/550; 
P_available_hp = Power*ones(size(V));

[LDmax, idx] = max(LD); 
V_LDmax = V(idx);
%% Stall speed (velocity where CL = CL_MAX)
V_stall = sqrt(2 * W / (rho * Sw * CL_MAX));

%% Maximum speed (velocity where Power Required <= Power Available)
V_max_idx = find(Power_Required_hp <= P_available_hp, 1, 'last');
if ~isempty(V_max_idx)
    V_max = V(V_max_idx);
else
    V_max = NaN;
end

%% Rate of climb (ft/s)
ROC = (P_available_hp - Power_Required_hp) * 550 ./ W;  % ft*lbf/s / lb = ft/s

% Rate of climb at stall
if ~isnan(V_stall)
    ROC_stall = ROC(find(V >= V_stall, 1, 'first'));
else
    ROC_stall = NaN;
end

%% ========================================================================
% 9) OUTPUT STORAGE
outputs = struct;
outputs.V = V;
outputs.CL = CL;
outputs.alpha = alpha;
outputs.alpha_t = alpha_t;
outputs.CLw = CLw;
outputs.CLt = CLt;
outputs.CDp = CDp;
outputs.CDi_w = CDi_w;
outputs.CDi_t = CDi_t;
outputs.CD_total = CD_total;
outputs.Dp = Dp;
outputs.Di = Di;
outputs.D = D;
outputs.LD = LD;
outputs.V_LDmax = V_LDmax;
outputs.LDmax = LDmax;
outputs.Power_Required_hp = Power_Required_hp;
outputs.P_available_hp = P_available_hp;
outputs.Re_w = Re_w;
outputs.Re_t = Re_t;
outputs.Re_f = Re_f;
outputs.lt_str = lt_str;
outputs.lt_tot = lt_tot;
outputs.lt_0 = lt_0;
outputs.VH_str = VH_str;
outputs.VH_tot = VH_tot;
outputs.VH_0 = VH_0;
outputs.hn = hn;
outputs.hcg_str = hcg_str;
outputs.hcg_tot = hcg_tot;
outputs.hcg_0 = hcg_0;
outputs.SM_str = SM_str;
outputs.SM_tot = SM_tot;
outputs.SM_0 = SM_0;
outputs.V_stall = V_stall;
outputs.V_max   = V_max;
outputs.ROC     = ROC;
outputs.ROC_stall = ROC_stall;   


if makePlots
    % Power
    figure; hold on
    plot(V, Power_Required_hp,'b-','LineWidth',2)
    plot(V, P_available_hp,'r--','LineWidth',2)
    xlabel('Velocity [ft/s]'); ylabel('Power [hp]')
    title('Power Required vs Velocity'); legend('Power Req','Power Avail'); grid on

    % Drag
    figure; hold on
    plot(V,Dp,'LineWidth',2); plot(V,Di,'LineWidth',2); plot(V,D,'LineWidth',2)
    plot(V_LDmax,D(idx),'ko','MarkerSize',8,'LineWidth',2)
    xlabel('Velocity [ft/s]'); ylabel('Drag [lb]')
    title('Drag Forces vs Velocity'); legend('Dp','Di','D','(L/D)_{max}'); grid on

    % L/D
    figure; hold on
    plot(V,LD,'LineWidth',2); plot(V_LDmax,LDmax,'ko','MarkerSize',8,'LineWidth',2)
    xlabel('Velocity [ft/s]'); ylabel('L/D'); title('L/D vs Velocity'); legend('L/D','(L/D)_{max}'); grid on
end

if makeTable
    alpha_deg = rad2deg(alpha); alpha_t_deg = rad2deg(alpha_t);
    T = table(V(:), CL(:), alpha_deg(:), alpha_t_deg(:), CLw(:), CLt(:), ...
        CDp(:), CDi_w(:), CDi_t(:), CD_total(:), Dp(:), Di(:), D(:), LD(:), ...
        Re_w(:), Re_t(:), Re_f(:), ...
        'VariableNames', {'V_ft_s','CL','alpha_deg','alpha_t_deg','CLw','CLt','CDp','CDi_w','CDi_t','CD_total','Dp_lb','Di_lb','D_lb','L_over_D','Re_w','Re_t','Re_f'});
    disp(T)
end

if makePrint_tail_Volume
    fprintf("\n==================== TAIL VOLUME ====================\n");
    fprintf("Tail AC station x_tac = %.4f ft\n", x_tac);
    fprintf("lt_str = %.4f ft --> VH_str = %.4f\n", lt_str, VH_str);
    fprintf("lt_tot = %.4f ft --> VH_tot = %.4f\n", lt_tot, VH_tot);
    fprintf("lt_0   = %.4f ft --> VH_0   = %.4f\n", lt_0,   VH_0);
end    

if makePrint_stability
    fprintf("\n==================== STABILITY ====================\n");
    fprintf("Neutral point hn = %.4f (fraction of MAC from wing LE)\n", hn);
    fprintf("hcg_str = %.4f  --> SM_str = %.4f (%.1f%% MAC)\n", hcg_str, SM_str, 100*SM_str);
    fprintf("hcg_tot = %.4f  --> SM_tot = %.4f (%.1f%% MAC)\n", hcg_tot, SM_tot, 100*SM_tot);
    fprintf("hcg_0   = %.4f  --> SM_0   = %.4f (%.1f%% MAC)\n", hcg_0, SM_0, 100*SM_0);
end
