%% Study of COSMIC monochromator
%%
% Computing the efficiency of gratings for COSMIC-U

%%
% requires Blazr.m
% May 2019, last updated August 2019
% awojdyla@lbl.gov

%% Grating parameters

% As of August 19, 2019
p_m = 24.301;
q_m = 7.573;

% pick a grating 
i_n = 1;

if i_n == 1
    lambda0_m = 1239/379.1*1e-9;
    k0 = 178.96e3 ;
    c = 1.632;
    material = 'Au';
elseif i_n == 2
    lambda0_m = 1239/806*1e-9;
    k0 = 287.44e3 ;
    c = 1.5772;
    material = 'Au';
    
elseif i_n == 3
    lambda0_m = 1239/1714.4*1e-9;
    k0 = 352.45e3 ;
    c = 1.7313;
    material = 'Rh';
    
end

%% trajectories

Es_eV = linspace(250,2500,100);
lambdas_m = 1.2398e-06./Es_eV;

% angles for all other wavelentgh
[alphas_rad, betas_rad] = Blazr.trajectory_vls(lambdas_m, lambda0_m, k0, p_m, q_m, c);

% compute the c-values and the demagnification
thetas_rad = (alphas_rad-betas_rad)/2;  % half-included angle
cs = cos(betas_rad)./cos(alphas_rad);   % c-values
magc = q_m./(cs*p_m);                    % demagnification

plot(Es_eV, alphas_rad*180/pi, Es_eV, abs(betas_rad)*180/pi, Es_eV, thetas_rad*180/pi,'k')
legend('entrance angle \alpha', 'exit angle |\beta|', 'half included angle \theta','location','Southeast')
xlabel('photon energy [eV]')
ylabel('angle [deg]')
title(sprintf('COSMIC-U trajectories (density=%1.1fl/mm, p=%1.1fm, q=%1.1fm, c=%1.2f)', k0*1e-3, p_m, q_m, c))
minb = round(min(-betas_rad*180/pi)*100)/100;
maxb = round(max(-betas_rad*180/pi)*100)/100;
mint = round(min(thetas_rad*180/pi)*100)/100;
maxt = round(max(thetas_rad*180/pi)*100)/100;
set(gca,'yTick',sort([88:0.5:90,minb,mint]))
xlim([250,2500])
grid on

fprintf('\nko=%1.1fl/mm\nmin theta = %1.2fdeg,\nmax theta = %1.2fdeg,\nmin beta = %1.2fdeg,\nmax beta = %1.2fdeg\n',...
    k0*1e-3, 90-mint, 90-maxt, 90-minb, 90-maxb)

%% Blazed grating efficiency
% (may take a while)

pitch_m = 1/k0;

% analytical blaze angle
theta_blaze_rad = (alpha0_rad + beta0_rad)/2;

% corresponding blazed grating height
thickness_m = 1/k0*tan(theta_blaze_rad);

etas_blaze = zeros(1,length(lambdas_m));
for i_l=1:length(lambdas_m)
    etas_blaze(i_l) = Blazr.efficiency_blazed(pitch_m, thickness_m, ...
        lambdas_m(i_l), pi/2-alphas_rad(i_l), material);
end

plot(Es_eV, etas_blaze)
xlabel('photon energy [eV]')
ylabel('efficiency')
title(sprintf('blazed grating; density=%1.1fl/mm, %s, thickness=%1.1fnm, c=%1.2f', k0*1e-3, material, thickness_m*1e9,c))
ylim([0 1])
grid on