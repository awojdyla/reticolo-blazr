%   1  D   exemple2_1D
%  VARIATION DE LA TRANSMISSION  FONCTION DE L'EPAISSEUR  

%clear;
LD=6;% longueur d'onde
D=10;% pas du reseau

teta0=0;nh=1;beta0=nh*sin(teta0*pi/180);
pol=-1; % 1:TE   -1:TM
parm=res0(pol);%  % initialisation des parametres par defaut
parm.sym.x=0;% utilisation de la symetrie


nn=10;% ordres de fourier 

% description des textures
textures{1}= 1;   
textures{2}= 1.5  ; 
textures{3}={[-4,4],[1.5,1]  };
% initialisation
aa=res1(LD,D,textures,nn,beta0,parm);
% definition du profil et calcul de la diffraction (on ne fait ce calcul qu'une fois)
figure;
t=[];epaisseur=[];
for h=linspace(0,40,50);
profil={[0,h,0] ,[1,3,2]  };
ef=res2(aa,profil);
% efficacite de transmission d'un incident du haut 
t=[t,ef.inc_top_transmitted.efficiency{0}];
epaisseur=[epaisseur,h];
plot(epaisseur,t);xlabel('epaisseur');ylabel('transmission dans l''ordre 0');title('VARIATION DE LA TRANSMISSION DANS L''ORDRE 0 FONCTION DE L''EPAISSEUR');pause(eps);
end

retio;



