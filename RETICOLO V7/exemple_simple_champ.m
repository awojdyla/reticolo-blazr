%%%%%%%%%%%%%%%%%%%%%%%
% EXEMPLE SIMPLE 2D   %
%%%%%%%%%%%%%%%%%%%%%%%

wavelength=8;% longueur d'onde
period=[10,15];% pas du reseau (meme unite que LD
n_incident_medium=1;%indice en haut
n_transmited_medium=1.5;%indice en haut
% incidence k_parallel(peut etre calcule par l'angle_theta dans le milieu incident) et delta0
angle_theta=10;k_parallel=n_incident_medium*sin(angle_theta*pi/180);
angle_delta=-20;

parm=res0;         % parametres par defaut
parm.res1.champ=1; % modification de la valeur par defaut pour le calcul du champ
nn=[3,2];  % ordres de fourier veut dire -nn(1) a nn(1) en x -nn(2) a nn(2) en y
% description des textures y compris le substrat et le superstrat et les milieux homogenes
texture=cell(1,3);
textures{1}= n_incident_medium;                   % milieu homogene 
textures{2}= n_transmited_medium;                 % milieu homogene  
textures{3}={n_incident_medium,[0,0,5,2,n_transmited_medium,1] }; % indice de base       inclusion 

aa=res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);

profile={[4.1,5.2,4.1],[1,3,2]};
deux_D=res2(aa,profile)

eff_TETM=search([-1,1],deux_D.TEinc_top_reflected,'efficiency')% <--efficacite (TE+TM) dans l'ordre -1,1 pour un incident TE du haut
eff_TE=search([-1,1],deux_D.TEinc_bottom_transmitted,'efficiency_TE') % <--efficacite TE dans l'ordre -1,1 pour un incident TE du bas

% calcul et trace des champs (coupe y=constante)  
x=linspace(-period(1)/2,period(1)/2,51);y=0;% x,y coordonnes du maillage en x et y (les cordonnes en z sont determinees par le calcul)
einc=[0,1];   %  composantes du champ E incident dans le repere u v  par defaut le champ incident vient du haut
parm.res3.trace=1; %trace automatique
parm.res3.npts=[50,50,50];
[e,z,o]=res3(x,y,aa,profile,einc,parm);

retio % efface les eventuels fichiers temporaires
 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXEMPLE SIMPLE CONIQUE %
%%%%%%%%%%%%%%%%%%%%%%%%%%

wavelength=8;% longueur d'onde
period=10;% pas du reseau (meme unite que LD
n_incident_medium=1;%indice en haut
n_transmited_medium=1.5;%indice en haut
% incidence k_parallel(peut etre calcule par l'angle angle_theta dans le milieu incident) et delta0
angle_theta0=10;k_parallel=n_incident_medium*sin(angle_theta0*pi/180);
angle_delta=-20;

parm=res0;         % parametres par defaut
parm.res1.champ=1; % modification de la valeur par defaut pour le calcul du champ
nn=5;  % ordres de fourier veut dire -nn a nn
% description des textures y compris le substrat et le superstrat et les milieux homogenes
texture=cell(1,3);
textures{1}= n_incident_medium;                   % milieu homogene 
textures{2}= n_transmited_medium;                   % milieu homogene  
textures{3}={[-2.5,2.5],[n_incident_medium,n_transmited_medium] }; % points de discontinuite indices a gauche 

aa=res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);

profile={[4.1,5.2,4.1],[1,3,2]};

conique=res2(aa,profile)

eff_TETM=search(1,conique.TEinc_top_reflected,'efficiency')% <--efficacite (TE+TM) dans l'ordre 1 pour un incident TE du haut
eff_TE=search(1,conique.TEinc_bottom_transmitted,'efficiency_TE') % <--efficacite TE dans l'ordre 1 pour un incident TE du bas


% calcul et trace des champs 
x=linspace(-period(1)/2,period(1)/2,51);% x coordonnes du maillage en x  (les cordonnes en z sont determinees par le calcul)
einc=[0,1];  %  composantes du champ E incident dans le repere u v  par defaut le champ incident vient du haut
parm.res3.trace=1; %trace automatique
parm.res3.npts=[50,50,50];
[e,z,o]=res3(x,aa,profile,einc,parm);
retio % efface les eventuels fichiers temporaires 
 
%%%%%%%%%%%%%%%%%%%%%%%%%
% EXEMPLE SIMPLE 1D TE  %
%%%%%%%%%%%%%%%%%%%%%%%%%

wavelength=8;% longueur d'onde
period=10;% pas du reseau (meme unite que LD
n_incident_medium=1;%indice en haut
n_transmited_medium=1.5;%indice en haut
% incidence k_parallel(peut etre calcule par l'angle angle_theta dans le milieu incident) et delta0
angle_theta0=-10;k_parallel=n_incident_medium*sin(angle_theta0*pi/180);

parm=res0(1);      % parametres par defaut et choix le la polarisation TE
parm.res1.champ=1; % modification de la valeur par defaut pour le calcul du champ
nn=40;  % ordres de fourier veut dire -nn a nn
% description des textures y compris le substrat et le superstrat et les milieux homogenes
texture=cell(1,3);
textures{1}= n_incident_medium;                   % milieu homogene 
textures{2}= n_transmited_medium;                   % milieu homogene  
textures{3}={[-2.5,2.5],[n_incident_medium,n_transmited_medium] }; % points de discontinuite indices a gauche 

aa=res1(wavelength,period,textures,nn,k_parallel,parm);

profile={[4.1,5.2,4.1],[1,3,2]};

un_D_TE=res2(aa,profile)
eff=search(1,un_D_TE.TEinc_top_reflected,'efficiency')% <--efficacite  dans l'ordre 1 pour un incident TE du haut

% calcul et trace des champs 
x=linspace(-period(1)/2,period(1)/2,51);% x coordonnes du maillage en x  (les cordonnes en z sont determinees par le calcul)
einc=1;  %  composantes du champ E incident sur uTE
parm.res3.trace=1; %trace automatique
parm.res3.npts=[50,50,50];
[e,z,o]=res3(x,aa,profile,einc,parm);

retio % efface les eventuels fichiers temporaires 

 
%%%%%%%%%%%%%%%%%%%%%%%%%
% EXEMPLE SIMPLE 1D TM  %
%%%%%%%%%%%%%%%%%%%%%%%%%

wavelength=8;% longueur d'onde
period=10;% pas du reseau (meme unite que LD
n_incident_medium=1;%indice en haut
n_transmited_medium=1.5;%indice en haut
% incidence k_parallel(peut etre calcule par l'angle angle_theta dans le milieu incident) et delta0
angle_theta0=-10;k_parallel=n_incident_medium*sin(angle_theta0*pi/180);

parm=res0(-1);      % parametres par defaut et choix le la polarisation TM
parm.res1.champ=1; % modification de la valeur par defaut pour le calcul du champ
nn=40;  % ordres de fourier veut dire -nn a nn
% description des textures y compris le substrat et le superstrat et les milieux homogenes
texture=cell(1,3);
textures{1}= n_incident_medium;                   % milieu homogene 
textures{2}= n_transmited_medium;                   % milieu homogene  
textures{3}={[-2.5,2.5],[n_incident_medium,n_transmited_medium] }; % points de discontinuite indices a gauche 

aa=res1(wavelength,period,textures,nn,k_parallel,parm);

profile={[4.1,5.2,4.1],[1,3,2]};
un_D_TM=res2(aa,profile)
eff=search(1,un_D_TM.TMinc_top_reflected,'efficiency')% <--efficacite  dans l'ordre 1 pour un incident TM du haut


% calcul et trace des champs 
x=linspace(-period(1)/2,period(1)/2,51);% x coordonnes du maillage en x  (les cordonnes en z sont determinees par le calcul)
einc=1;  %  composantes du champ H incident sur uTE
parm.res3.trace=1; %trace automatique
parm.res3.npts=[50,50,50];
[e,z,o]=res3(x,aa,profile,einc,parm);


retio % efface les eventuels fichiers temporaires 

