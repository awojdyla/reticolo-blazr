%   2  D    exemple_general_conique
%  EXEMPLE GENERAL CONIQUE


%  toutes les longueurs doivent etre exprimees dans la meme unite que la longueur d'onde 
%                   les angles sont en degres 
% simplification du programme 2D au cas conique
clear;


% parametres generaux
%........................
LD=10;% longueur d'onde
D=17;% pas du reseau


% incidence ro(peut etre calcule par l'angle teta dans le milieu incident) et delta0
teta0=10;nh=1;ro=nh*sin(teta0*pi/180);
delta0=20;
delta0=0;



% parametre technique
%........................
nn=5;  % ordres de fourier


textures=cell(1,5);
% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}={ 1};   
textures{2}={ 1.5  };
textures{3}={ [-4,4] ,  [1,1.5]    };   
textures{4}={[-1,-.5,.5,1],[2,1.3,1.5,1.3]  };   

%textures{5}={ inf,[0,4,1]           };%  metal infiniment conducteur

textures{5}={ [-1,1],        [1.5,2.5]            };
%            profil      indices agauche


%*****************************************************************
% ETAPE D'INITIALISATION  calcul du descripteur des textures
%******************************************************************
%  cette etape doit etre refaite chaque fois que l'on change LD D nn les symetries ro delta0
%  par contre elle ne depend pas du profil du reseau construit a partir des textures (sequences hauteurs) 

% symetrie

parm=res0;

parm.sym.x=0; % x=x0 :plan de symetrie seulement si beta0(1)=0    par defaut parm.sym.x=[]  pas de symetrie en x 
%parm.sym.pol=-1;% 1 TE  -1:TM

% parametres facultatifs controlant cette initialisation
parm.res1.trace=1; % trace des textures
%  parm.res1.xlimite=[0,D(1)]; parm.res1.nx=500; % parametres pour le trace des textures
%  parm.res1.calcul=0;% pas calcul si on en veut que verifier les textures
 parm.res1.champ=1; % si on veut un calcul soigne du champ
%  parm.res1.ftemp=1;  % utilisation de fichiers temporaires
%  parm.res1.fperm='aaa'; % utilisation de fichiers permanents
%  parm.res1.sog=0; % matrices G


%cal=1;if cal==1;
aa=res1(LD,D,textures,nn,ro,delta0,parm);%  aa contient toutes les informations pour la suite
%save test aa;%else;load test aa;end;


%***********************************************************************************************************
% ETUDE DE RESAUX
%***********************************************************************************************************


% description du profil du reseau
%************************************
%        hauteurs=[0,1,.5,.6,0];           sequence=[1,3,2,4,2];      profil={hauteurs,sequence};
%              hauteurs            numeros des textures  correspondantes
%    la premiere et la derniere  texture est une texture hogogene eventuellement de hau teur nulle(superstrat et substrat)
%      ou        profil={{hauteurs,sequence}, nfois} si on veut repeter nfois le motif
% on donne ici 4 exemples

profil0={[0,0], [1,2]  };
profil1={      [0,1,.5,.6,0]    [1,4,5,4,2]       };
profil2={ {0,1} , {  [1,.5,.6] , [3,2,4] ,2 },   {[2,0],[5,2]}    };
profil3={ {1,1} , {  [1,.5,.6] ,[3,2,4]  ,2 },   {[2,1],[5,2]}    }; % idem a 2 avec des hauteurs en haut et en bas
profil4={ {0,1} , {  [1,.5,.6] ,[3,2,4] ,2 },   {[6,3] ,[5,2] } ,{0,2} };


% verification des profils par trace de coupes
%**********************************************

x=linspace(-D(1)/2,D(1)/2,50);% coordonnes du maillage en x  z est determinee par le programme
parm=res0;
parm.res3.cale=[] ;%signifie que l'on ne calcule pas le champ
parm.res3.trace=1 ;%trace automatique


[tab1,z,o]=res3(x,aa,profil1,[1,1],parm); % profil1 
[tab1,z,o]=res3(x,aa,profil2,[1,1],parm); % profil2 
[tab1,z,o]=res3(x,aa,profil3,[1,1],parm); % profil3 
[tab4,z,o]=res3(x,aa,profil4,[1,1],parm); % profil4 

% on peut aussi se contenter d'imprimer les tableaux tab qui decrivent les profils
parm.res3.caltab=1; % on ne calcule que tab
[tab1,z,o]=res3(x,aa,profil1,[1,1],parm);tab1, % profil1 
[tab2,z,o]=res3(x,aa,profil2,[1,1],parm);tab2, % profil1 
[tab3,z,o]=res3(x,aa,profil3,[1,1],parm);tab3, % profil1 

%*******************************************************************
% etude de la diffraction d'une onde plane incidente par un reseau
%*******************************************************************

ef=res2(aa,profil1);


% la structure ef est la 'carte d'identitee des ordres incidents et diffractes  propagatifs 
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
x=linspace(-D(1)/2,D(1)/2,50);;% coordonnes en x  (les cordonnes en z sont determinees par le calcul)
einc=[1,1];   %  composantes du champ e incident dans le repere u v  par defaut le champ incident vient du haut
parm=res0;parm.res3.trace=1;
[e,z,o]=res3(x,aa,profil3,einc,parm);

% calcul des champs       ex(z,x)=e(:,:,1),ey(z,x)=e(:,:,2),ez(z,x)=e(:,:,3)
%                         hx(z,x)=e(:,:,4),hy(z,x)=e(:,:,5),hz(z,x)=e(:,:,6)
%  de l'indice complexe du reseau     n(z,x)=o
%  et des valeurs de z  le bas de la structure etudiee est z=0


retio % efface les eventuels fichiers temporaires 