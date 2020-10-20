%   1  D     exemple_1D_metal
%  EXEMPLE METAL PARFAIT


retio;
% parametres generaux
%........................
LD=2*pi;;% longueur d'onde
D=5.5;% pas du reseau
% incidence beta0(peut etre calcule par l'angle teta dans le milieu incident) 
teta0=0;nh=1;beta0=nh*sin(teta0*pi/180);
nn=5;  % ordres de fourier
% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= inf;   
textures{2}={ inf,[0,1,1]     };
textures{3}={ inf,[0,5,1]     };
textures{4}= 1.5 ;
textures{5}= 1;
textures{6}={ [-D/4,D/4],[1,1.5]     };
textures{7}={inf, [D/2,D/2,1]     };
%*****************************************************************
% ETAPE D'INITIALISATION  calcul du descripteur des textures
%******************************************************************
pol=1;% 1:TE   -1:TM
parm=res0(pol);
parm.sym.x=0; % x=x0 :plan de symetrie seulement si beta0=0    par defaut parm.sym.x=[]  pas de symetrie  
parm.res1.trace=1; % trace des textures  
parm.res1.champ=0; % si on ne veut pas calculer les champs
aa=res1(LD,D,textures,nn,beta0,parm);%  aa contient toutes les informations pour la suite

% description du profil du reseau
profil={      [1.5,.1,3.27,1]  ,  [1,2,3,4]       };
profil1={      [1,.5,D,1] ,   [5,7,6,1]       };
profil2={      [1,D,.5,1]  ,  [1,6,7,5]       };

%*******************************************************************
% calcul et trace des champs
%*******************************************************************
x=linspace(-D,D,101);% coordonnes en x  (les cordonnees en z sont determinees par le calcul)

parm=res0(pol);  % initialisation des parametres par defaut
parm.res3.sens=-1;
parm.res3.npts=[10,5,50,5];
parm.res3.trace=1 ;%trace automatique
[e,z,o]=res3(x,aa,profil,1,parm);

parm.res3.npts=[10,10,80,5];
parm.res3.sens=1;
[e,z,o]=res3(x,aa,profil1,1,parm);


parm.res3.npts=[5,80,10,10];
parm.res3.sens=-1;
[e,z,o]=res3(x,aa,profil2,1,parm);


retio % efface les eventuels fichiers temporaires 
















