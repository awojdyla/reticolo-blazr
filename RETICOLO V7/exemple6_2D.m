%   2  D    exemple6_2D
%  VARIATION DE LA REFLEXION  FONCTION DE TETA  

clear;retio
figure;
r=[];teta=[];

for teta0=sort([linspace(0,20,21),linspace(0,89,21)]); 
LD=10;
D=[10,10];% pas du reseau

nh=1;ro=nh*sin(teta0*pi/180);
delta0=0;
parm=res0;
parm.sym.y=0;parm.sym.x=0;parm.sym.pol=1;% utilisation des symetries si l'incidence le permet

nn=[5,5];% ordres de fourier 

% description des textures
textures{1}={ 1};   
textures{2}={ 1.5  }; 
textures{3}={ 1.5   ,  [0,0,5,8,  1,  1]     };
% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);
% definition du profil et calcul de la diffraction (on doit refaire ce calcul pour chaque LD)

profil={[0,20,0]  ,[1,3,2]  };
    
ef=res2(aa,profil);
% TE incident du haut,reflexion dans l'ordre 0
r=[r,ef.TEinc_top_reflected.efficiency_TE{0}];teta=[teta,teta0];
plot(teta,r);xlabel('teta degres');ylabel('reflexion');title('VARIATION DE LA REFLEXION FONCTION DE TETA');pause(eps);
end

retio;



