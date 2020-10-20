%   2  D    exemple_general_2D
% EXEMPLE GENERAL

%  toutes les longueurs doivent etre exprimees dans la meme unite que la longueur d'onde 
%                   les angles sont en degres 


clear;retio;

% parametres generaux
%........................
LD=8;% longueur d'onde
D=[10,15];% pas du reseau

% incidence ro(peut etre calcule par l'angle teta dans le milieu incident) et delta0
teta0=10;nh=1;ro=nh*sin(teta0*pi/180);
delta0=-20;

% parametre technique
%........................
nn=[3,2];  % ordres de fourier veut dire -nn(1) a nn(1) en x -nn(2) a nn(2) en y


% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= 1;   
textures{2}= 1.5; 
textures{3}={ 1   ,  [0,0,5,2,  1.5,  1] ,  [0,0,1,10,  1.5,  1]          };   
textures{4}={ 2.5      [0,0,5,5,  1.1,  1]      };   
textures{5}={ 2.5      [0,0,8,6,  1,  2]                                                              };   
textures{6}={ 1.5  ,      [0,0,6,6,      1,      5],             [0,0,2,2,  1.5,  1]                          };
%         indice de base       inclusion 1                          inclusion 2
%                         ellipse approchee par 5 rectangles      rectangle 
%                         de centre 0,0 de grands axes 6 6      de centre 0,0 de cotes 2 2
%                               et d' indice 1                        et d'indice 1.5
%      pour les profils 1D en x invariant en y   
textures{7}={[-.5,1,5,3],[2,1.3,1.5,2.1]   };
textures{8}={[-.5,1,5,3],[2,1.3,1.5,2.1] , [0,0,6,6,  1,  5] ,  [0,0,2,2,  1.5,  1]    };
textures{9}={[-.5,1,5,3].',[2,1.3,1.5,2.1].' , [0,0,6,6,  1,  5] ,  [0,0,2,2,  1.5,  1]    };
%                  profil      indices

%*****************************************************************
% ETAPE D'INITIALISATION  calcul du descripteur des textures
%******************************************************************
%  cette etape doit etre refaite chaque fois que l'on change LD D nn les symetries ro delta0
%  par contre elle ne depend pas du profil du reseau construit a partir des textures (sequences hauteurs) 

% parametres facultatifs controlant cette initialisation
parm=res0;  % initialisation des parametres par defaut

parm.sym.x=0; % x=x0 :plan de symetrie seulement si beta0(1)=0    par defaut parm.sym.x=[]  pas de symetrie en x 
parm.sym.y=0; % y=y0 :plan de symetrie seulement si beta0(2)=0   par defaut parm.sym.y=[]  pas de symetrie en y
parm.sym.pol=-1;% 1 TE  -1:TM choix de la polarisation si la symetrie est prise en compte





  parm.res1.trace=1; % trace des textures
% parm.res1.xlimite=[-D(1),D(1)]; parm.res1.nx=500; % parametres pour le trace des textures en x
%  parm.res1.ylimite=[-D(2),D(2)]; parm.res1.nx=500; % parametres pour le trace des textures en y
%  parm.res1.calcul=0;% pas calcul si on en veut que verifier les textures
  parm.res1.champ=1; % si on veut un calcul soigne du champ
%  parm.res1.fperm='aaa'; % utilisation de fichiers permanents
%  parm.res1.sog=0; % matrices G

%cal=0;if cal==1;
[aa,nef]=res1(LD,D,textures,nn,ro,delta0,parm);
%  aa contient toutes les informations pour la suite
%  nef{ii} est le vecteur des indices efficaces de textures{ii}
%save test aa;else;load test aa;end;

%***********************************************************************************************************
% ETUDE DE RESAUX
%***********************************************************************************************************


% description du profil du reseau
%************************************
%   sequence=[1,3,2,4,2];        hauteurs=[0,1,.5,.6,0];       profil={sequence,hauteurs};
%   numeros des textures   hauteurs correspondantes
%    la premiere et la derniere  texture est une texture homogene eventuellement de hauteur nulle(superstrat et substrat)
%      ou        profil={{sequence,hauteurs}, nfois} si on veut repeter nfois le motif
% on donne ici 4 exemples

profil0={[0,0],[1,2]  };
profil1={    [0,1,.5,.6,0] , [1,4,6,4,2]           };
profil2={ {0,1} , {  [1,.5,.6]  , [3,2,4] ,2 },   {[2,0] ,[6,2]}    };
profil3={ {1,1} , {  [1,.5,.6]  , [3,2,4],2 },   {[2,1]  ,[6,2]}    }; % idem a 2 avec des hauteurs en haut et en bas
profil4={ {0,1} ,  {  [1,.5,.6]  , [3,2,4] ,2 },   {[1,2]  ,[6,3]}  ,{  [1,0,.2]  , [5,2,3] ,2 },{0,2} };


% verification des profils par trace de coupes
%**********************************************

x=linspace(-D(1)/2,D(1)/2,51);y=0;% coordonnes du maillage en x y ( z est determinee par le programme )
parm=res0;
parm.res3.cale=[] ;%signifie que l'on ne calcule pas le champ
parm.res3.trace=1 ;%trace automatique


[tab1,z,o]=res3(x,y,aa,profil1,[1,1],parm); % profil1 
[tab2,z,o]=res3(x,y,aa,profil2,[1,1],parm); % profil2 
[tab3,z,o]=res3(x,y,aa,profil3,[1,1],parm); % profil3 
% on peut aussi se contenter de calculer les tableaux tab qui decrivent les profils
parm.res3.caltab=1; % on ne calcule que tab
[tab1,z,o]=res3(x,y,aa,profil1,[1,1],parm);tab1, % profil1 
[tab2,z,o]=res3(x,y,aa,profil2,[1,1],parm);tab2, % profil1 
[tab3,z,o]=res3(x,y,aa,profil3,[1,1],parm);tab3, % profil1 

%*******************************************************************
% etude de la diffraction d'une onde plane incidente par un reseau
%*******************************************************************

ef=res2(aa,profil1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% la structure ef est la 'carte d'identitee'  des ondes planes incidentes et diffractees  propagatives 
%  et fournit en outre les amplitudes diffractees par une onde incidente 
%
%     certains tableaux peuvent etre vides  (notamment a cause des symetries)
%    et certains champs inexistants( cas 1 D)
%
% 
%                           TEinc_top: [1x1 struct]
%                 TEinc_top_reflected: [1x1 struct]
%               TEinc_top_transmitted: [1x1 struct]
% 
%                        TEinc_bottom: [1x1 struct]
%              TEinc_bottom_reflected: [1x1 struct]
% ef         TEinc_bottom_transmitted: [1x1 struct]
% 
%                           TMinc_top: [1x1 struct]
%                 TMinc_top_reflected: [1x1 struct]
%               TMinc_top_transmitted: [1x1 struct]
% 
%                        TMinc_bottom: [1x1 struct]
%              TMinc_bottom_reflected: [1x1 struct]
%            TMinc_bottom_transmitted: [1x1 struct]
%  
% 
% 
%                                  
% avec par exemple:
%                          theta: [5x1 double]
%                          delta: [5x1 double]
%                              K: [5x3 double] vecteur d'onde unitaire
%
%                     efficiency: [5x1 double]
%                  efficiency_TE: [5x1 double]
%                  efficiency_TM: [5x1 double]
%                   amplitude_TE: [5x1 double]
%                   amplitude_TM: [5x1 double]
%
%  ef.TEinc_top_transmitted    E: [5x3 double]  composantes de E a l'origine 
%                              H: [5x3 double]  composantes de H a l'origine 
%
%                 PlaneWave_TE_E: [5x3 double] composantes de l'onde plane TE  (TM)
%                 PlaneWave_TE_H: [5x3 double] dont la composante du vecteur de 
%                                              poynting selon z est 1/2
%
%                PlaneWave_TE_Eu: [5x2 double]  memes ondes 
%                PlaneWave_TE_Hu: [5x2 double]  dans le repere u_TM ,u_TE
%
%                 PlaneWave_TM_E: [5x3 double]
%                 PlaneWave_TM_H: [5x3 double]
%
%                PlaneWave_TM_Eu: [5x2 double]
%                PlaneWave_TM_Hu: [5x2 double]



enerh=sum(ef.TEinc_top_reflected.efficiency);enerb=sum(ef.TEinc_top_transmitted.efficiency);
disp(['TE incident du haut : energie en haut  ',num2str(enerh),' energie en bas ',num2str(enerb),' pertes ',num2str(1-enerh-enerb)])

enerh=sum(ef.TEinc_bottom_reflected.efficiency);enerb=sum(ef.TEinc_bottom_transmitted.efficiency);
disp(['TE incident du bas : energie en haut  ',num2str(enerh),' energie en bas ',num2str(enerb),' pertes ',num2str(1-enerh-enerb)])


enerh=sum(ef.TMinc_top_reflected.efficiency);enerb=sum(ef.TMinc_top_transmitted.efficiency);
disp(['TM incident du haut : energie en haut  ',num2str(enerh),' energie en bas ',num2str(enerb),' pertes ',num2str(1-enerh-enerb)])

enerh=sum(ef.TMinc_bottom_reflected.efficiency);enerb=sum(ef.TMinc_bottom_transmitted.efficiency);
disp(['TM incident du bas : energie en haut  ',num2str(enerh),' energie en bas ',num2str(enerb),' pertes ',num2str(1-enerh-enerb)])



%*******************************************************************
% calcul et trace des champs
%*******************************************************************


% coupe y=constante  
x=linspace(-D(1)/2,D(1)/2,51);y=0;% x,y coordonnes du maillage en x et y (les cordonnes en z sont determinees par le calcul)
einc=[1,1];   %  composantes du champ e incident dans le repere u v  par defaut le champ incident vient du haut
parm=res0; %parametres par defaut
parm.res3.trace=1 ;%trace automatique

[e,z,o]=res3(x,y,aa,profil3,einc,parm);

% calcul des champs       ex(z,x,y)=e(:,:,:,1),ey(z,x,y)=e(:,:,:,2),ez(z,x,y)=e(:,:,:,3)
%                         hx(z,x,y)=e(:,:,:,4),hy(z,x,y)=e(:,:,:,5),hz(z,x,y)=e(:,:,:,6)
%  de l'indice complexe du reseau     n(z,x,y)=o(:,:,:)
%  et des valeurs de z  le bas de la structure etudiee est z=0
 
% coupe x=constante
x=0;y=linspace(-D(2)/2,D(2)/2,51);einc=[1,1];
parm=res0; %parametres par defaut
parm.res3.npts=5;% nombre de points par tranche
parm.res3.sens=0; % champ incident de dessous
parm.res3.trace=1 ;%trace automatique
parm.res3.champs=[1,2,i];
[e,z,o]=res3(x,y,aa,profil3,einc,parm); 



% on peut par exemple calculer le champ sur un plan a une hauteur z0 dans une tranche

x=linspace(-D(1)/2,D(1)/2,51);y=linspace(-D(2)/2,D(2)/2,51);einc=[1,1];
z0=.32;
profil={  [0,1,2-z0,0,z0,1,0] , [1,4,3,3,3,4,2] };einc=[1,1];
parm=res0; %parametres par defaut
parm.res3.npts=[0,0,0,1,0,0,0];% nombre de points par tranche
parm.res3.trace=1;parm.res3.champs=[i,1,2];%trace automatique de Ex Hx avec contour de l'objet

[e,z,o]=res3(x,y,aa,profil,einc,parm); 


retio % efface les eventuels fichiers temporaires 