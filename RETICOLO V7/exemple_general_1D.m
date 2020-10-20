%   1  D    exemple_general_1D
%  EXEMPLE GENERAL



%  toutes les longueurs doivent etre exprimees dans la meme unite que la longueur d'onde 
%                   les angles sont en degres 
%  programme 1D 


% parametres generaux
%........................
LD=10;% longueur d'onde
D=7;% pas du reseau


% incidence beta0(peut etre calcule par l'angle teta dans le milieu incident) 
teta0=10;nh=1;beta0=nh*sin(teta0*pi/180);



% parametre technique
%........................
nn=5;  % ordres de fourier



% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= 1;   

textures{2}= 1.5  ; 
textures{3}={ [-4,4] ,  [1,1.5]    };   
textures{4}={[-1,-.5,.5,1],[2,1.3,1.5,1.3]  };   
textures{5}={ 2.5                           };

textures{6}={ [-1,1],        [1.5,2.5]      };
%textures{6}={ inf,[3.5,3,1.5]     };
%  TEXTURES{ }= { inf, [cx1,dx1,ni1],[cx2,dx2,ni2],..} pour le metal 

%            profil      indices a gauche


%*****************************************************************
% ETAPE D'INITIALISATION  calcul du descripteur des textures
%******************************************************************
%  cette etape doit etre refaite chaque fois que l'on change LD D nn les symetries beta0
%  par contre elle ne depend pas du profil du reseau construit a partir des textures (sequences hauteurs) 

pol=1;% 1:TE   -1:TM
parm=res0(pol);  % initialisation des parametres par defaut

parm.res1.trace=1;  
%parm.res1.champ=0;% si on ne veut pas calculer les champs  
parm.sym.x=0; % x=x0 :plan de symetrie seulement si beta0=0    par defaut parm.sym.x=[]  pas de symetrie  


aa=res1(LD,D,textures,nn,beta0,parm);%  aa contient toutes les informations pour la suite


%***********************************************************************************************************
% ETUDE DE RESEAUX
%***********************************************************************************************************


% description du profil du reseau
%************************************
%        hauteurs=[0,1,.5,.6,0];           sequence=[1,3,2,4,2];      profil={hauteurs,sequence};
%              hauteurs            numeros des textures  correspondantes
%    la premiere et la derniere  texture est une texture hogogene eventuellement de hauteur nulle(superstrat et substrat)
%      ou        profil={{hauteurs,sequence}, nfois} si on veut repeter nfois le motif
% on donne ici 4 exemples

profil0={[0,0], [1,2]  };
profil1={      [0,1,.5,.6,0]    [1,4,6,4,2]       };
profil2={ {0,1} , {  [1,.5,.6] , [3,2,4] ,2 },   {[2,0],[6,2]}    };
profil3={ {1,1} , {  [1,.5,.6] ,[3,2,4]  ,2 },   {[2,1],[6,2]}    }; % idem a 2 avec des hauteurs en haut et en bas
profil4={ {0,1} ,  {  [1,.5,.6] , [3,2,4] ,2 },   {[6,3] ,[1,2] }  ,{  [1,0,.2] , [5,2,3] ,2 },{0,2} };
ef1=res2(aa,profil1,parm);

% verification des profils par trace de coupes
%**********************************************

x=linspace(0,D(1),50);% coordonnes du maillage en x  z est determinee par le programme
parm=res0(pol);%parm.res2.result=0;  % initialisation des parametres par defaut

parm.res3.cale=[] ;%signifie que l'on ne calcule pas le champ
parm.res3.trace=1 ;%trace automatique

[tab1,z,o]=res3(x,aa,profil1,[1,1],parm); % profil1 
[tab1,z,o]=res3(x,aa,profil2,[1,1],parm); % profil2 
[tab1,z,o]=res3(x,aa,profil3,[1,1],parm); % profil3
[tab4,z,o]=res3(x,aa,profil4,[1,1],parm); % profil4 

% on peut aussi se contenter de calculer les tableaux tab qui decrivent les profils
parm.res3.caltab=1; % on ne calcule que tab
[tab1,z,o]=res3(x,aa,profil1,[1,1],parm);tab1, % profil1 
[tab2,z,o]=res3(x,aa,profil2,[1,1],parm);tab2, % profil2 
[tab3,z,o]=res3(x,aa,profil3,[1,1],parm);tab3, % profil2

%*******************************************************************
% etude de la diffraction d'une onde plane incidente par un reseau
%*******************************************************************

ef=res2(aa,profil1,parm);

% la structure ef est la 'carte d'identitee des ordres incidents et diffractes  propagatifs 
% 
%                           inc_top: [1x1 struct]
%                inc_top_reflected: [1x1 struct]
%               inc_top_transmitted: [1x1 struct]
% ef
%                       inc_bottom: [1x1 struct]
%              inc_bottom_reflected: [1x1 struct]
%           inc_bottom_transmitted: [1x1 struct]
% 
%  
% 
% 
%                                  
% avec par exemple:
%                          theta: [5x1 double]
%                              K: [5x3 double] vecteur d'onde unitaire
%
%                     efficiency: [5x1 double]
%                     amplitude: [5x1 double]
%
%  ef.inc_top_transmitted      E: [5x3 double]  composantes de E a l'origine 
%                              H: [5x3 double]  composantes de H a l'origine 
%
%                    PlaneWave_E: [5x3 double] composantes de l'onde plane TE  (TM)
%                    PlaneWave_H: [5x3 double] dont la composante du vecteur de 
%                                              poynting selon z est 1/2


                        

disp('incident du haut')
disp(rettexte('en haut:',ef.inc_top_reflected.efficiency));
disp(rettexte('en bas:',ef.inc_top_transmitted.efficiency))
disp(rettexte('total:',sum(ef.inc_top_reflected.efficiency)+sum(ef.inc_top_transmitted.efficiency)))

% amplitudes diffractees par l' incident du bas
disp('incident du bas')
disp(rettexte('en haut:',ef.inc_bottom_reflected.efficiency));
disp(rettexte('en bas:',ef.inc_bottom_transmitted.efficiency))
disp(rettexte('total:',sum(ef.inc_bottom_reflected.efficiency)+sum(ef.inc_bottom_transmitted.efficiency)))

%*******************************************************************
% calcul et trace des champs
%*******************************************************************


% coupe  
x=linspace(-D(1)/2,D(1)/2,51);% coordonnes en x  (les cordonnees en z sont determinees par le calcul)

parm=res0(pol);
parm.res3.sens=1;
parm.res3.trace=1 ;%trace automatique

[e,z,o]=res3(x,aa,profil3,1,parm);

% calcul des champs    en TE   Ey(z,x,y)=e(:,:,1), Hx(z,x)=e(:,:,2), Hz(z,x,y)=e(:,:,3)
%                      en TM   Hy(z,x,y)=e(:,:,1), Ex(z,x)=e(:,:,2),Ez(z,x,y)=e(:,:,3)
%  de l'indice complexe du reseau     n(z,x)=o
%  et des valeurs de z  (le bas de la structure etudiee est z=0 )
%  calcul du champ    par defaut le champ incident vient du haut d'amplitude 1

% et pour l'incident du bas:
parm=res0(pol);
parm.res3.sens=-1;
parm.res3.trace=1;parm.res3.champs=[i,1,2,3];%trace automatique avec coutours de l'objet
parm.res3.champs=[i,1,2,3];

[e,z,o]=res3(x,aa,profil3,1,parm);

retio % efface les eventuels fichiers temporaires 