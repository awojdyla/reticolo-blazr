%   1  D   exemple9_1D
% METAL INFINIMENT CONDUCTEUR RECHERCHE  D'UNE RESONANCE, TRACE DU CHAMP A LA RESONANCE

clear;
LD=15;% longueur d'onde
D=10;% pas du reseau

nn=10 ;% ordres de fourier 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= 1 ;   
textures{2}= inf  ; 
textures{3}={ inf   ,  [0,3,  1]     };
textures{4}={ inf   ,  [0,8,  1]     };

% initialisation
parm=res0(1); % initialisation des parametres par defautpolarisation TE
parm.sym.x=0;;% utilisation de la symetrie

aa=res1(LD,D,textures,nn,0,parm);


% recherche de la resonance en h comme pole de h-->R
er=inf;iter=0;z0=20;z=[];zz=[];
while(er>eps)&(iter<30);
profil={[0,1,z0,0] ,[1,3,4,2]  };
ef=res2(aa,profil);
zz0=1/ef.inc_top_reflected.amplitude{0};
iter=iter+1;z=[z;z0];zz=[zz;zz0];[z0,z,zz,er]=retcadillac(z,zz);
end;

h=real(z0);


rr=[];
hh=h+linspace(-.1,.1,30);
for hhh=hh;
profil={[0,1,hhh,0] ,[1,3,4,2]  };
ef=res2(aa,profil);
rr=[rr,ef.inc_top_reflected.amplitude{0}];
end;
figure;plot(hh,real(rr),hh,imag(rr));xlabel('hauteur');legend('real R','imag R');title(['resonance pour h=',num2str(h)]);pause(eps);

% coupe y=0  
x=linspace(-1.5*D(1),1.5*D(1),150);

parm.res3.npts=[10,5,40,10];% nombre de points par tranche
parm.res3.trace=1 ;%trace automatique


profil={[5,1,h,5] , [1,3,4,2]  };
res3(x,aa,profil,1,parm);


retio;