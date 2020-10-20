%   1  D    exemple6_1D
%  VARIATION DE LA TRANSMISSION DANS L'ORDRE 1 FONCTION DE TETA  

clear;retio

LD=6;% longueur d'onde
D=10;% pas du reseau

pol=-1; % 1:TE   -1:TM


nn=10;% ordres de fourier 

% description des textures
textures{1}= 1;   
textures{2}= 1.5  ; 
textures{3}={[-4,4],[1.5,1]  };





figure;
t=[];teta=[];

for teta0=linspace(-89,89,101); 
nh=1;beta0=nh*sin(teta0*pi/180);

 % initialisation
parm=res0(pol);  % initialisation des parametres par defaut
parm.sym.x=0;% utilisation de la symetrie si l'incidence le permet
aa=res1(LD,D,textures,nn,beta0,parm);
    
% definition du profil et calcul de la diffraction (on doit refaire ce calcul pour chaque LD)
profil={[0,20,0]  ,[1,3,2]  };
ef=res2(aa,profil);

t=[t,ef.inc_top_transmitted.efficiency{1}];teta=[teta,teta0];
plot(teta,t);xlabel('teta degres');ylabel('reflexion');title('VARIATION DE LA TRANSMISSION DANS L''ORDRE 1 FONCTION DE TETA');pause(eps);
end

retio;



