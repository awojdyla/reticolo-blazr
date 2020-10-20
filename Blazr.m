classdef Blazr < handle
    % BLAZR Grating monochromator tool
    
    % The software is provided "as is",
    % without warranty of any kind, express or implied
    % A Wojdyla, ALS/ALS-U, LBNL
    % awojdyla@lbl.gov
    % May 2019
    
    % add encircled energy, strehl, propagation, fwhm
    
    % Download RETICOLO v7 or v8
    %
    % (previous) https://www.lp2n.institutoptique.fr/Membres-Services/Responsables-d-equipe/LALANNE-Philippe
    % https://www.lp2n.institutoptique.fr/light-complex-nanostructures
    
    properties (Constant)
        c   = 299792458         % speed of light [m/s]
        h   = 6.626070040*1e-34 % Planck constanct [J.s]
        eV  = 1.6021766208*1e-19% elementaty charge x 1V [J.s]
        Ag  = 6.0221408571e23   % Avogadro number [mol-1]
        deg = 180/pi            % degree/rad conversion
        sigma_fhwm = 1/(2*sqrt(2*log(2))); % gaussian fwhm to variance
        in_m = 25.4e-3          % inch to meter conversion
    end
    
    % autoload script
    %     if exist('MIP','class')~=8
    %         websave('MIP.m','https://raw.githubusercontent.com/awojdyla/mip/master/MIP.m')
    %     end
    
    methods(Static)
        
        % beam size
       
        function [alphas_rad, betas_rad] = trajectory_vls(lambdas_m, lambda0_m, k0, p_m, q_m, c)
            % TRAJECTORY_VLS Calculate the monochromator trajectory for a VLS
            %   Blazr.trajectory_vls(lambdas_m, lambda0_m, k0, p_m, q_m, c)
            %
            % See also BLAZR.TRAJECTORY_C
            
            % diffracted order
            m = 1;
            
            % Energies
            Es_eV = 1.2398e-06./lambdas_m;
            
            beta0 = @(alpha_rad) asin(m*lambda0_m*k0-sin(alpha_rad));
            fun0 = @(alpha_rad) sin(alpha_rad) + sin(-acos(c*cos(alpha_rad))) - m*lambda0_m*k0;
            
            % angles at optimized wavelength
            alpha0_rad = fzero(fun0, pi/2);
            beta0_rad = beta0(alpha0_rad);
            
            % linear coefficient
            b2 = 1/(2*m*lambda0_m*k0)*(cos(alpha0_rad).^2/p_m+cos(beta0_rad).^2/q_m);
            
            alphas_rad = zeros(1, length(Es_eV));
            betas_rad  = zeros(1, length(Es_eV));
            cs         = zeros(1, length(Es_eV));
            for i_l = 1:length(Es_eV)
                lambda_m = lambdas_m(i_l);
                beta = @(alpha_rad) asin(m*lambda_m*k0-sin(alpha_rad));
                fun2 = @(alpha_rad) cos(alpha_rad).^2/p_m+cos(beta(alpha_rad)).^2/q_m - m*2*b2*lambda_m*k0;
                alphas_rad(i_l) = fzero(fun2, pi/2);
                betas_rad( i_l)  = beta(alphas_rad(i_l));
                cs(i_l) = cos(betas_rad(i_l))/cos(alphas_rad(i_l));
            end
            
        end
        
        
        function [b2, b3, b4] = VLS_parameters(lambda_m, k0, alpha_rad, beta_rad, p_m, q_m)
            % VLS_PARAMETERS Gives the VLS parameters for a specific setting
            %   [b2, b3, b4] = VLS_parameters(lambda_m, k0, alpha_rad, beta_rad, p_m, q_m)
            %
            %   based on doi.org/10.1016/j.nima.2004.09.007
            %
            % See also ...
            
            warning('be careful, not fully tested!')
            
            % linear coefficient
            b2 = 1/(2*m*lambda_m*k0)*(cos(alpha_rad).^2/p_m+cos(beta_rad).^2/q_m);
            b3 = 1/(2*m*lambda_m*k0)*(cos(alpha_rad).^2.*sin(alpha_rad)/p_m...
                -cos(beta_rad).^2.*sin(beta_rad)/q_m);
            b4 = 1/(8*m*lambda_m*k0)*( cos(alpha_rad).^2.*(4*sin(alpha_rad).^2-cos(alpha_rad).^2)/p_m ...
                +cos(beta_rad).^2.*(4*sin(beta_rad).^2-cos(beta_rad).^2)/q_m);
            
        end
        
        
        function eta = efficiency_lamellar(pitch_m, ...
                thickness_m, wavelength_m, grazing_angle_rad, material)
            % EFFICIENCY_LAMELLAR Diffraction efficiency of a lamellar grating.
            %
            %   eta = Blazr.efficiency_lamellar(pitch_m, ...
            %            thickness_m, wavelength_m, grazing_angle_rad, material)
            %
            % See also Blazr.efficiency_blazed
            
            % requires reticolo and list of material indices
            
            
            % fixed parameters
            
            % polarisation
            pol  = 1; % 1:TE   -1:TM
            % number of Fourier orders
            nn = 10;
            
            if strcmp(material,'Rh')
                load('n_rh.mat','E_eV','n_c')
            elseif strcmp(material,'Au')
                load('n_gold.mat','E_eV','n_c')
            elseif strcmp(material,'Ni')
                load('n_ni.mat','E_eV','n_c')
            else
                load('n_gold.mat','E_eV','n_c')
                warning('wrong material name, using Gold')
            end
            
            addpath(genpath('RETICOLO V7'))
            retio
            
            % grating pitch
            p_nm  = pitch_m*1e9;
            
            % angle of incidence
            theta0_deg=(pi/2-grazing_angle_rad)*180/pi;
            % thickness
            th_nm = thickness_m*1e9;
            
            beta0 = sin(theta0_deg*pi/180);
            
            parm = res0(pol);  % parameter init
            
            % wavelength
            lambda_nm = wavelength_m*1e9;
            
            % complex index of gold
            %n_au = 1-2.1062E-03+1i*1.0306E-03;
            n_au = conj(interp1(E_eV, n_c, 1239/lambda_nm));
            
            % description of textures
            textures{1} = 1; % top medium
            textures{2} = n_au; % bottom medium
            textures{3} ={[-p_nm/2,0],[1, n_au]};
            
            % initialisation
            parm.res1.trace= 0;
            aa = res1(lambda_nm, p_nm, textures, nn, beta0, parm);
            
            profile = {[0,th_nm,0] ,[1,3,2]};
            
            ef = res2(aa, profile);
            % grating efficiency
            eta = ef.inc_top_reflected.efficiency{-1};
            
            retio
            
        end
        
        function eta = efficiency_blazed(pitch_m, ...
                thickness_m, wavelength_m, grazing_angle_rad, material)
            % EFFICIENCY_BLAZED Diffraction efficiency of a blazed grating.
            %
            %   eta = Blazr.efficiency_blazed(pitch_m, ...
            %            thickness_m, wavelength_m, grazing_angle_rad, material)
            %
            % See also Blazr.efficiency_blazed
            
            % clear;
            addpath(genpath('RETICOLO V7'))
            retio
            
            
            % number of Fourier orders
            nn = 10;
            % material
            if strcmp(material,'Rh')
                load('n_rh.mat','E_eV','n_c')
            elseif strcmp(material,'Au')
                load('n_gold.mat','E_eV','n_c')
            elseif strcmp(material,'Ni')
                load('n_ni.mat','E_eV','n_c')
            else
                load('n_gold.mat','E_eV','n_c')
                warning('wrong material name, using Gold')
            end
            
            % number of blazed layers
            N_blaze_layers = 20;
            
            % grating pitch
            p_nm  = pitch_m*1e9;
            
            % angle of incidence
            theta0_deg = (pi/2-grazing_angle_rad)*180/pi;
            % thickness
            th_nm = thickness_m*1e9;
            
            %beta0 = th_nm*sin(theta0_deg*pi/180);
            beta0 = sin(theta0_deg*pi/180);
            
            % polarisation
            pol = -1; % 1:TE   -1:TM
            parm = res0(pol);  % parameter init
            
            % wavelength
            lambda_nm = wavelength_m*1e9;
            
            % complex index of gold
            n_mat = conj(interp1(E_eV, n_c, 1239/lambda_nm));
            %n_mat = 1-2.1062E-03+1i*1.0306E-03;
            
            textures = cell(1,N_blaze_layers+2);
            
            textures{1} = 1; % top medium
            textures{N_blaze_layers+2} = n_mat; % bottom medium
            
            A1 = zeros(1,N_blaze_layers);
            A2 = zeros(1,N_blaze_layers);
            for i_l = 1:N_blaze_layers
                textures{i_l+1} ={[-(i_l-1)*p_nm/N_blaze_layers,0],[1, n_mat]};
                A1(i_l) = th_nm/(N_blaze_layers);
                A2(i_l) = i_l+1;
            end
            
            profile = {[0,A1,0] ,[1,A2,(N_blaze_layers+2)]};
            
            % initialisation
            parm.res1.trace= 0;
            aa = res1(lambda_nm, p_nm, textures, nn, beta0, parm);
            
            % show if needed
            %x=linspace(-p_nm,p_nm,501);% on trace 2 periodes
            %parm.res3.trace=1 ; % trace automatique
            %parm.res3.cale=[];
            %parm.res3.npts=[10,80,10];
            %[e,z,o]=res3(x,aa,profile,1,parm);
            %axis square
            %set(gcf,'WindowStyle','docked')
            
            ef = res2(aa, profile, parm);
            % grating efficiency
            eta = ef.inc_top_reflected.efficiency{-1};
        end
        
        function eta = efficiency_lamellar_duty(pitch_m, ...
                thickness_m, wavelength_m, grazing_angle_rad, material, duty)
            % EFFICIENCY_LAMELLAR Diffraction efficiency of a lamellar grating.
            %
            %   eta = Blazr.efficiency_lamellar(pitch_m, ...
            %            thickness_m, wavelength_m, grazing_angle_rad, material)
            %
            % See also Blazr.efficiency_blazed
            
            % requires reticolo and list of material indices
            
            
            % fixed parameters
            
            % polarisation
            pol  = 1; % 1:TE   -1:TM
            % number of Fourier orders
            nn = 10;
            
            if strcmp(material,'Rh')
                load('n_rh.mat','E_eV','n_c')
            elseif strcmp(material,'Au')
                load('n_gold.mat','E_eV','n_c')
            elseif strcmp(material,'Ni')
                load('n_ni.mat','E_eV','n_c')
            else
                load('n_gold.mat','E_eV','n_c')
                warning('wrong material name, using Gold')
            end
            
            addpath(genpath('RETICOLO V7'))
            retio
            
            % grating pitch
            p_nm  = pitch_m*1e9;
            
            % angle of incidence
            theta0_deg=(pi/2-grazing_angle_rad)*180/pi;
            % thickness
            th_nm = thickness_m*1e9;
            
            beta0 = sin(theta0_deg*pi/180);
            
            parm = res0(pol);  % parameter init
            
            % wavelength
            lambda_nm = wavelength_m*1e9;
            
            % complex index of gold
            %n_au = 1-2.1062E-03+1i*1.0306E-03;
            n_au = conj(interp1(E_eV, n_c, 1239/lambda_nm));
            
            % description of textures
            textures{1} = 1; % top medium
            textures{2} = n_au; % bottom medium
            textures{3} ={[-p_nm*duty,0],[1, n_au]};
            
            % initialisation
            parm.res1.trace= 0;
            aa = res1(lambda_nm, p_nm, textures, nn, beta0, parm);
            
            profile = {[0,th_nm,0] ,[1,3,2]};
            
            ef = res2(aa, profile);
            % grating efficiency
            eta = ef.inc_top_reflected.efficiency{-1};
            
            retio
            
        end
        
        
        function eta = efficiency_blazed_higher(pitch_m, ...
                thickness_m, wavelength_m, grazing_angle_rad, material, order)
            % EFFICIENCY_BLAZED Diffraction efficiency of a blazed grating.
            %
            %   eta = Blazr.efficiency_blazed(pitch_m, ...
            %            thickness_m, wavelength_m, grazing_angle_rad, material)
            %
            % See also Blazr.efficiency_blazed
            
            % clear;
            addpath(genpath('RETICOLO V7'))
            retio
            
            
            % number of Fourier orders
            nn = 10;
            % material
            if strcmp(material,'Rh')
                load('n_rh.mat','E_eV','n_c')
            elseif strcmp(material,'Au')
                load('n_gold.mat','E_eV','n_c')
            elseif strcmp(material,'Ni')
                load('n_ni.mat','E_eV','n_c')
            else
                load('n_gold.mat','E_eV','n_c')
                warning('wrong material name, using Gold')
            end
            
            % number of blazed layers
            N_blaze_layers = 20;
            
            % grating pitch
            p_nm  = pitch_m*1e9;
            
            % angle of incidence
            theta0_deg = (pi/2-grazing_angle_rad)*180/pi;
            % thickness
            th_nm = thickness_m*1e9;
            
            %beta0 = th_nm*sin(theta0_deg*pi/180);
            beta0 = sin(theta0_deg*pi/180);
            
            % polarisation
            pol = -1; % 1:TE   -1:TM
            parm = res0(pol);  % parameter init
            
            % wavelength
            lambda_nm = wavelength_m*1e9;
            
            % complex index of gold
            n_mat = conj(interp1(E_eV, n_c, 1239/lambda_nm));
            %n_mat = 1-2.1062E-03+1i*1.0306E-03;
            
            textures = cell(1,N_blaze_layers+2);
            
            textures{1} = 1; % top medium
            textures{N_blaze_layers+2} = n_mat; % bottom medium
            
            A1 = zeros(1,N_blaze_layers);
            A2 = zeros(1,N_blaze_layers);
            for i_l = 1:N_blaze_layers
                textures{i_l+1} ={[-(i_l-1)*p_nm/N_blaze_layers,0],[1, n_mat]};
                A1(i_l) = th_nm/(N_blaze_layers);
                A2(i_l) = i_l+1;
            end
            
            profile = {[0,A1,0] ,[1,A2,(N_blaze_layers+2)]};
            
            % initialisation
            parm.res1.trace= 0;
            aa = res1(lambda_nm, p_nm, textures, nn, beta0, parm);
            
            % show if needed
            %x=linspace(-p_nm,p_nm,501);% on trace 2 periodes
            %parm.res3.trace=1 ; % trace automatique
            %parm.res3.cale=[];
            %parm.res3.npts=[10,80,10];
            %[e,z,o]=res3(x,aa,profile,1,parm);
            %axis square
            %set(gcf,'WindowStyle','docked')
            
            ef = res2(aa, profile, parm);
            % grating efficiency
            eta = ef.inc_top_reflected.efficiency{-order};
        end
    end
end