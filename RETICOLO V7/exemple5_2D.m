%   2  D    exemple5_2D
%  VARIATION DE LA TRANSMISSION DANS l'ORDRE (-1,0)  FONCTION DE LAMBDA 

clear;retio
figure;
t=[];lambda=[];

for LD=linspace(9,11,21);% longueur d'onde

D=[10,10];% pas du reseau

teta0=0;nh=1;ro=nh*sin(teta0*pi/180);
delta0=0;
parm=res0;
parm.sym.x=0;parm.sym.y=0;parm.sym.pol=1;% utilisation de 2 symetries

nn=[5,5];% ordres de fourier 

% description des textures
textures{1}= 1;   
textures{2}=1.5; 
textures{3}={ 1.5   ,  [0,0,5,8,  1,  1]     };
% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);
% definition du profil et calcul de la diffraction (on doit refaire ce calcul pour chaque LD)

profil={[0,20,0] ,[1,3,2]  };
    
ef=res2(aa,profil);

% transmission dans l''ordre (-1,0)
t=[t,ef.TEinc_top_transmitted.efficiency_TE{-1,0}];lambda=[lambda,LD];
plot(lambda,t);xlabel('longueur d''onde');ylabel('transmission dans l''ordre (-1,0)');title('VARIATION DE LA TRANSMISSION DANS l''ORDRE (-1,0) FONCTION DE LAMBDA');pause(eps);
end

retio;



