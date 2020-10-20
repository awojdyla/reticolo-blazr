%   1  D    exemple_1D_pertes
% CALCUL DES PERTES 1D

clear;
LD=8;6;% longueur d'onde
D=10;% pas du reseau

teta0=0;
nh=1.2;beta0=nh*sin(teta0*pi/180);
nb=1.5;
for pol=[-1,1]; % 1:TE   -1:TM
parm=res0(pol);
parm.sym.x=0;% utilisation de la symetrie
nn=50; % ordres de fourier 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}=nh;   
textures{2}=nb; 

for cas=1:3;
switch cas
case 1;textures{3}={[-1,1],[.1+5i,1]  };% dielectrique
case 2;textures{3}={inf, [-D/4,D/4,.1+5i] ,[D/4,D/4,.1+5i]   };% metal electrique
case 3;textures{3}={-inf, [0,D/2,.1+5i]   };% metal magnetique
end;

profil={[3,.3,1.4] ,[1,3,2]  };parm.res3.npts=[[0,10,0];[1,8,1]];
% initialisation
aa=res1(LD,D,textures,nn,beta0,parm);
ef=res2(aa,profil);

x=[-D/2,D/2];
parm.res3.trace=0 ; % trace automatique

parm.res3.gauss_x=20;% 10 par defaut

% [e,z,o]=res3(x,aa,profil,1,parm);
parm.res3.sens=1;% incident du haut
if pol==1;inc=ef.inc_top.PlaneWave_E(2);else;inc=ef.inc_top.PlaneWave_H(2);end;
[e,z,o,w,PP,P,p,XX]=res3(x,aa,profil,inc,parm);PP
bilan_energie=sum(ef.inc_top_reflected.efficiency)+sum(ef.inc_top_transmitted.efficiency)+sum(PP)/(.5*D)-1
figure;retcolor(XX,z,p);

parm.res3.sens=-1;% incident du bas
if pol==1;inc=ef.inc_bottom.PlaneWave_E(2);else;inc=ef.inc_bottom.PlaneWave_H(2);end;
[e,z,o,w,PP,P,p,XX,wx]=res3(x,aa,profil,inc,parm);PP
bilan_energie=sum(ef.inc_bottom_reflected.efficiency)+sum(ef.inc_bottom_transmitted.efficiency)+sum(PP)/(.5*D)-1
figure;pcolor(XX,z,p);shading flat;

if pol==1;
p1=pi/LD*imag(o.^2).*abs(e(:,:,1)).^2;
else;
p1=pi/LD*imag(o.^2).*sum(abs(e(:,:,2:3)).^2,3);
end;    
P1=(p1*wx(:)).';
retcompare(p,p1),retcompare(P,P1)
erreur_pertes=max(abs(p(:)-p1(:)))/(max(abs(p(:)))+max(abs(p1(:)))),max(abs(P-P1))/(max(abs(P))+max(abs(P1))),

retio;
end; % cas
end; % pol
