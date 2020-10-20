%   1  D    exemple3_1D
%  TRACE DU CHAMP

clear;
LD=6;% longueur d'onde
D=10;% pas du reseau

teta0=0;nh=1;beta0=nh*sin(teta0*pi/180);
pol=-1; % 1:TE   -1:TM
parm=res0(pol);
parm.sym.x=0;% utilisation de la symetrie
nn=20;% ordres de fourier 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}=1;   
textures{2}= 1.5; 
textures{3}={[-4,4],[1.5,1]  };
% initialisation
aa=res1(LD,D,textures,nn,beta0,parm);

profil={[3,20,3] ,[1,3,2]  };
%%
x=linspace(-D,D,501);% on trace 2 periodes
parm.res3.trace=1 ; % trace automatique
parm.res3.npts=[10,80,10];

[e,z,o]=res3(x,aa,profil,1,parm);


retio;