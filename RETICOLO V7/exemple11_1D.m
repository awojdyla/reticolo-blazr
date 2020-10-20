%   1  D    exemple11_1D
%  RESEAU 1D  TRANSMISSION EXCEPTIONNELLE 
%  transmission dans l'ordre 0 fonction de l'epaisseur
% (comparer a exemple11_2D)

clear;

t=[];h=[];

LD=6;% longueur d'onde
D=5;% pas du reseau
pol=-1;
parm=res0(pol);  % initialisation des parametres par defaut
parm.sym.x=0;% utilisation de la symetrie
pol=-1;

nn=10;% ordres de fourier 

% description des textures
nm=0.030465+3.2792i;
textures{1}= 1;   
textures{2}={[-1,1],[ 0.030465+3.2792i,1]  }; 

% initialisation
beta0=0;
aa=res1(LD,D,textures,nn,beta0,parm);
for epaisseur=linspace(0,LD/2,1000);
% definition du profil et calcul de la diffraction 
profil={[0,epaisseur,0] ,[1,2,1]  };

[ef,tab]=res2(aa,profil);
% transmission dans l'ordre 0 (incident du haut)
h=[h,epaisseur];
t=[t,ef.inc_top_transmitted.efficiency{0}];
end;
figure;
plot(h,t);xlabel('epaisseur');title('RESEAU 1D METALLIQUE TRANSMISSION EXCEPTIONNELLE');ylabel('efficacitee transmise');pause(eps);

retio