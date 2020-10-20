%% Comparison between Reticolo and PCGrate
%%
% Computing the efficiency of a blazed grating to compare with results
% calculated by Dima (at constant angle) using PC Grate

%%
% September 2019
% awojdyla@lbl.gov

%% Grating paramaters

p_m = 24.301;
q_m = 7.573;

lambda0_m = 1239/806*1e-9;
g = 287e3 ;
c = 1.577;
m = 1;
material = 'Au';

Es_eV = linspace(100,2500,240);
lambdas_m = 1.2398e-06./Es_eV;

% analytical blaze angle
%theta0_blaze_rad = (alpha0_rad + beta0_rad)/2;
%thickness_m = 1/g*tan(theta0_blaze_rad);

%% Trajectories

% for the design wavelength
x0= (-(2*m*lambda0_m*g)+sqrt((2*m*lambda0_m*g)^2-4*(c^2-1)*(1-m^2*lambda0_m^2*g^2-c^2)))/(2*(c^2-1));
y0 = -sqrt(1-c^2*(1-x0^2));
alpha0_rad = asin(x0);
beta0_rad = asin(y0);

% defocus coefficient
b2 = (cos(alpha0_rad).^2/p_m + cos(beta0_rad).^2/q_m)/(2*g*lambda0_m);

alphas_rad = zeros(1,length(Es_eV));
betas_rad  = zeros(1,length(Es_eV));
thetas_rad = zeros(1,length(Es_eV));

for i_e=1:length(Es_eV)
    
    lambda_m  = lambdas_m(i_e);

    p2 = -(1/p_m+1/q_m);
    p1 = 2*m*g*lambda_m/q_m;
    p0 = (1/p_m + 1/q_m - m^2*g^2*lambda_m^2/q_m - 2*m*g*b2*lambda_m);

    x = (-p1-sqrt(p1^2-4*p2*p0))/(2*p2);
    alpha_rad = asin(x);

    y = x-m*g*lambda_m;
    beta_rad = -asin(y);

    theta_rad = (alpha_rad - beta_rad)/2;

    alphas_rad(i_e) = alpha_rad;
    betas_rad(i_e)  = beta_rad;
    thetas_rad(i_e) = theta_rad;
end

% shorthand
%[alphas_rad, betas_rad] = Blazr.trajectory_vls(lambdas_m, lambda0_m, g, p_m, q_m, c);
%thetas_rad = (alphas_rad-betas_rad)/2;  % half-included angle

%alphap_rad = sqrt(2*1*lambda0_m/(1/g*(c^2-1)));

plot(Es_eV, alphas_rad*180/pi, Es_eV, abs(betas_rad)*180/pi, Es_eV, thetas_rad*180/pi,'k')
legend('entrance angle \alpha', 'exit angle |\beta|', 'half included angle \theta','location','Southeast')
xlabel('photon energy')
ylabel('angle [deg]')
title(sprintf('COSMIC-U trajectories (density=%1.1fl/mm, p=%1.1fm, q=%1.1fm, c=%1.2f)', g*1e-3, p_m, q_m, c))
%set(gca,'yTick',87.5:0.5:90)
grid on



%% Blazed grating efficiency at constant angle
% (may take a while)

i_0 = find(Es_eV>=806,1);
alpha0_rad = alphas_rad(i_0);
beta0_rad  =  betas_rad(i_0);

pitch_m = 1/g;

% analytical blaze angle
theta0_blaze_rad = (alpha0_rad + beta0_rad)/2;
thetaN_blaze_rad = 0.2*pi/180;
% corresponding blazed grating height
thickness_m = 1/g*tan(thetaN_blaze_rad);

etas_blaze2 = zeros(1,length(lambdas_m));
etas_blaze4 = zeros(1,length(lambdas_m));
etas_blaze5 = zeros(1,length(lambdas_m));
etas_blaze6 = zeros(1,length(lambdas_m));

fprintf('\nmight take a few minutes on a slow computer\n')
for i_l=1:length(lambdas_m)
    string_status = sprintf('calculating wavelengh %i/%i', i_l, length(lambdas_m));
    fprintf(string_status)
    
    etas_blaze2(i_l) = Blazr.efficiency_blazed(pitch_m, 1/g*tan(0.2*pi/180), ...
        lambdas_m(i_l), pi/2-alpha0_rad, material);
    etas_blaze4(i_l) = Blazr.efficiency_blazed(pitch_m, 1/g*tan(0.4*pi/180), ...
        lambdas_m(i_l), pi/2-alpha0_rad, material);
    etas_blaze5(i_l) = Blazr.efficiency_blazed(pitch_m, 1/g*tan(0.5*pi/180), ...
        lambdas_m(i_l), pi/2-alpha0_rad, material);
    etas_blaze6(i_l) = Blazr.efficiency_blazed(pitch_m, 1/g*tan(0.6*pi/180), ...
        lambdas_m(i_l), pi/2-alpha0_rad, material);
    % completion status
    
    n_char = length(string_status);
    for i_c=1:n_char
        fprintf('\b')
    end
end

plot(Es_eV, etas_blaze2, '.-c',...
     Es_eV, etas_blaze4, '.-b',...
     Es_eV, etas_blaze5, '.-g',...
     Es_eV, etas_blaze6, '.-r')
xlabel('photon energy [eV]')
ylabel('efficiency')
title(sprintf('RETICOLO: blazed grating; density=%1.1fl/mm, c=%1.2f, material:%s', g*1e-3, c, material))
legend('\theta_b=0.2deg','\theta_b=0.4deg','\theta_b=0.5deg','\theta_b=0.6deg')
ylim([0 0.5])
set(gca,'ytick',0:0.05:1)
set(gca,'xtick',200:200:2500)
grid on
drawnow

%% Blazed grating efficiency along the trajectory


i_0 = find(Es_eV>=806,1);
alpha0_rad = alphas_rad(i_0);
beta0_rad  =  betas_rad(i_0);

pitch_m = 1/g;

etas_blaze2 = zeros(1,length(lambdas_m));
etas_blaze4 = zeros(1,length(lambdas_m));
etas_blaze5 = zeros(1,length(lambdas_m));
etas_blaze6 = zeros(1,length(lambdas_m));
for i_l=1:length(lambdas_m)
    etas_blaze2(i_l) = Blazr.efficiency_blazed(pitch_m, 1/g*tan(0.2*pi/180), ...
        lambdas_m(i_l), pi/2-alphas_rad(i_l), material);
    etas_blaze4(i_l) = Blazr.efficiency_blazed(pitch_m, 1/g*tan(0.4*pi/180), ...
        lambdas_m(i_l), pi/2-alphas_rad(i_l), material);
    etas_blaze5(i_l) = Blazr.efficiency_blazed(pitch_m, 1/g*tan(0.5*pi/180), ...
        lambdas_m(i_l), pi/2-alphas_rad(i_l), material);
    etas_blaze6(i_l) = Blazr.efficiency_blazed(pitch_m, 1/g*tan(0.6*pi/180), ...
        lambdas_m(i_l), pi/2-alphas_rad(i_l), material);
end

plot(Es_eV, etas_blaze2, '.-c',...
     Es_eV, etas_blaze4, '.-b',...
     Es_eV, etas_blaze5, '.-g',...
     Es_eV, etas_blaze6, '.-r')
xlabel('photon energy [eV]')
ylabel('efficiency')
title(sprintf('RETICOLO: blazed grating; density=%1.1fl/mm, c=%1.2f, material:%s', g*1e-3, c, material))
ylim([0 0.5])
set(gca,'ytick',0:0.05:1)
set(gca,'xtick',200:200:2500)
grid on
legend('\theta_b=0.2deg','\theta_b=0.4deg','\theta_b=0.5deg','\theta_b=0.6deg')
