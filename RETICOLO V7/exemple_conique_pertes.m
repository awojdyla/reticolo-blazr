%   2  D    exemple_conique_pertes
% CALCUL DES PERTES CONIQUE

clear;
LD=8;% longueur d'onde
D=10;% pas du reseau

teta0=30;nh=1.5;ro=nh*sin(teta0*pi/180);
nb=1.2;

delta0=20;
parm=res0;

parm.res1.champ=1; % si on veut un calcul soigne du champ

nn=20;% ordres de fourier 

for cas=1:3;

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= nh;   
textures{2}= nb; 
switch cas;
case 1;textures{3}={ [-2.5,2.5] , [.1+5i,1]};% dielectriques
case 2;textures{3}={ inf,[ 0 , D/2, .1+5i]};% metaux electriques
case 3;textures{3}={ -inf,[ 0 , D/2, .1+5i]};% metaux magnetiques
end;
% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);


x=[-D(1)/2,D(1)/2];

parm.res3.npts=[[0,10,0];[1,5,1]];% nombre de points par tranche
profil={[5,.2,5] ,[1,3,2]  };
[ef,tab]=res2(aa,profil);

parm.res3.trace=0;
einc=ef.TEinc_top.PlaneWave_TE_Eu;[e,z,o,w,PP,P,p,XX,wx]=res3(x,aa,profil,einc,parm);  %   incident  TE
bilan_energie=sum(ef.TEinc_top_reflected.efficiency)+sum(ef.TEinc_top_transmitted.efficiency)+sum(PP)/(.5*D)-1
figure;pcolor(XX,z,p);shading flat;
p1=pi/LD*imag(o.^2).*sum(abs(e(:,:,1:3)).^2,3);
P1=(reshape(p1,length(z),[])*wx(:)).';
erreur_pertes=max(abs(p(:)-p1(:)))/(max(abs(p(:)))+max(abs(p1(:)))),max(abs(P-P1))/(max(abs(P))+max(abs(P1))),

einc=ef.TMinc_top.PlaneWave_TM_Eu;[e,z,o,w,PP,P,p,XX,wx]=res3(x,aa,profil,einc,parm);  %   incident  TE
bilan_energie=sum(ef.TMinc_top_reflected.efficiency)+sum(ef.TMinc_top_transmitted.efficiency)+sum(PP)/(.5*D)-1
figure;pcolor(XX,z,p);shading flat;
p1=pi/LD*imag(o.^2).*sum(abs(e(:,:,1:3)).^2,3);
P1=(reshape(p1,length(z),[])*wx(:)).';
erreur_pertes=max(abs(p(:)-p1(:)))/(max(abs(p(:)))+max(abs(p1(:)))),max(abs(P-P1))/(max(abs(P))+max(abs(P1))),

retio;

end;

