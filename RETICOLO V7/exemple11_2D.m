%   2  D    exemple11_2D
%  RESEAU 1D traite come conique avec delta=0
%  transmission dans l'ordre 0 fonction de l'epaisseur
%    (comparer a exemple11_1D)
clear;

t=[];h=[];

LD=6;% longueur d'onde
D=5;% pas du reseau
parm=res0;
parm.sym.y=0;%  symetries par rapport a y=0  
parm.sym.x=0;%  symetries par rapport a x=0  
parm.sym.pol=-1;  % 1 TE  -1 TM

nn=10;% ordres de fourier 

% description des textures
nm=0.030465+3.2792i;
textures{1}= 1;   
textures{2}={[-1,1],[ 0.030465+3.2792i,1]  }; 

% initialisation
ro=0;delta0=0;
aa=res1(LD,D,textures,nn,ro,delta0,parm);
for epaisseur=linspace(0,LD/2,1000);
% definition du profil et calcul de la diffraction 
profil={[0,epaisseur,0] ,[1,2,1]  };

[ef,tab]=res2(aa,profil);
h=[h,epaisseur];
t=[t,ef.TMinc_top_transmitted.efficiency_TM{0}];% TM incident du haut,transmission TM
end;
figure;
plot(h,t);xlabel('epaisseur');title('RESEAU 1D METALLIQUE TRANSMISSION EXCEPTIONNELLE');ylabel('efficacitee transmise');pause(eps);

retio