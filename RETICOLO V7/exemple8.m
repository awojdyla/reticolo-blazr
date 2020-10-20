% exemple8
%  UTILISATION DE FICHIERS PERMANENTS

clear;
LD=6;% longueur d'onde
D=[10,10];% pas du reseau
teta0=0;nh=1;ro=nh*sin(teta0*pi/180);
delta0=0;
parm=res0;parm.res2.result=0;
parm.sym.x=0;parm.sym.y=0;parm.sym.pol=1;% utilisation de 2 symetries
parm.res1.fperm='fich'; %  le resultat est mis sur un fichier permanent  de nom 'fich123..'qu'il ne faut pas oublier d'effacer ensuite
  
nn=[5,5];% ordres de fourier 

% description des textures
textures{1}= 1;   
textures{2}= 1.5  ; 
textures{3}={ 1.5   ,  [0,0,5,8,  1,  1]     };
% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);

aa=res1(LD,D,textures,nn,ro,delta0,parm);

clear;retio;
aa='fich';
parm=res0;
profil={[0,20,0] ,[1,3,2] };
ef=res2(aa,profil);
ef
retio;



