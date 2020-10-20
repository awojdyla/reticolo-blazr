%% Blazr Mimimal example

% grating pitch
pitch_m = 5.5878e-06;
% grating material  (choose from 'Au', 'Ni', 'Rh')
material = 'Au';

% wavelenth (380eV)
lambda0_m = 3.2683e-09;
% incidence angle
alpha0_rad = 1.5443;
% exit angle (+1 order)
beta0_rad = -1.5275;

%% efficiency for a lamellar grating

thickness_lam_m = 50e-9;
eta_lamellar = Blazr.efficiency_lamellar(pitch_m, thickness_lam_m, ...
        lambda0_m, pi/2-alpha0_rad, material);

fprintf('efficiency for a lamellar grating is %02.1f percent\n', eta_lamellar*100)
    
%% efficiency for a blazed grating
    
% analytical blaze angle
theta_blaze_rad = (alpha0_rad + beta0_rad)/2;
% corresponding blazed grating height
thickness_blaze_m = pitch_m*tan(theta_blaze_rad);

eta_blazed = Blazr.efficiency_blazed(pitch_m, thickness_blaze_m, ...
        lambda0_m, pi/2-alpha0_rad, material);
    
fprintf('efficiency for a   blazed grating is %02.1f percent\n', eta_blazed*100)
