%   2  D    exemple7_2D
%  VARIATION DE LA REFLEXION  FONCTION DE DELTA 

clear;retio
figure;
r=[];delta=[];
teta0=60;
nh=1;ro=nh*sin(teta0*pi/180);
LD=10;
D=[10,10];% pas du reseau
% ici on ne peut pas utiliser les symetries
nn=[2,2];% ordres de fourier 

% description des textures
textures{1}= 1;   
textures{2}= 1.5 ; 
textures{3}={ 1.5   ,  [0,0,5,8,  1,  1]     };
profil={[0,20,0],[1,3,2]  };
for delta0=linspace(0,180,101); 
aa=res1(LD,D,textures,nn,ro,delta0);
% initialisation
parm=res0;
ef=res2(aa,profil,parm);

% onde TE incidente du haut,reflexion TE dans l'ordre 0 
r=[r,ef.TEinc_top_reflected.efficiency_TE{0}];delta=[delta,delta0];
plot(delta,abs(r).^2);xlabel('delta degres');ylabel('reflexion');title('VARIATION DE LA REFLEXION FONCTION DE DELTA');pause(eps);
end
retio;



