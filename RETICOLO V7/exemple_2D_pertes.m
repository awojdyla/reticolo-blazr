%   2  D    exemple_2D_pertes
% CALCUL DES PERTES 2D


clear;
LD=8;% longueur d'onde
D=[10,8];% pas du reseau

teta0=10;nh=1.5;ro=nh*sin(teta0*pi/180);
nb=1.2;
delta0=0;
parm=res0;
parm.sym.x=0;parm.sym.y=0;parm.sym.pol=1;% utilisation de 2 symetries
parm.res1.champ=1; % si on veut un calcul soigne du champ

nn=[5,5];% ordres de fourier 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}=nh;   
textures{2}= nb; 
for cas=2;1:3;
switch cas
case 1;textures{3}={ .1+5i   ,  [0,0,9,1,  1,  1]     }; % dielectrique
case 2;textures{3}={ inf   ,  [-D(1)/4,0,D(1)/4,D(2)/2,  .1+5i] ,  [D(1)/4,0,D(1)/4,D(2)/2,  .1+5i]     }; % metal electrique
case 3;textures{3}={-inf   ,  [0, 0, D(1)/4,D(2)/2,  .1+5i]     }; % metal magnetique
end;

% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);

profil={[1,1.5,1] ,[1,3,2]  };
ef=res2(aa,profil);    

% coupe y=0  
x=linspace(0,D(1),150);y=0;
if parm.sym.pol==1
einc=ef.TEinc_top.PlaneWave_TE_Eu;   %  composantes du champ e incident  TE
else;
einc=ef.TMinc_top.PlaneWave_TM_Eu;   %  composantes du champ e incident  TE
end   
    
    
parm.res3.npts=[[0,10,0];[1,2,1]];% nombre de points par tranche methode de Gauss
% 
x=[-D(1)/2,D(1)/2];y=[-D(2)/2,D(2)/2];parm.res3.trace=0;
[e,z,o,w,PP,P,p,XX,YY,wxy]=res3(x,y,aa,profil,einc,parm);
[prv,iY0]=min(abs(YY));figure;pcolor(XX,z,p(:,:,iY0));shading flat;title('pertes');
if parm.sym.pol==1;
bilan_energie=sum(ef.TEinc_top_reflected.efficiency)+sum(ef.TEinc_top_transmitted.efficiency)+sum(PP)/(.5*prod(D))-1
else;
bilan_energie=sum(ef.TMinc_top_reflected.efficiency)+sum(ef.TMinc_top_transmitted.efficiency)+sum(PP)/(.5*prod(D))-1
end 


p1=pi/LD*imag(o.^2).*sum(abs(e(:,:,:,1:3)).^2,4);
P1=(reshape(p1,length(z),[])*wxy(:)).';
erreur_pertes=max(abs(p(:)-p1(:)))/(max(abs(p(:)))+max(abs(p1(:)))),max(abs(P-P1))/(max(abs(P))+max(abs(P1))),

retio;
end;
