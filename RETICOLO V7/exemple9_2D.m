%   2  D    exemple9_2D
% METAL INFINIMENT CONDUCTEUR RECHERCHE D' UNE RESONANCE, TRACE DU CHAMP A LA RESONANCE 

clear;
LD=15;% longueur d'onde
D=[10,10];% pas du reseau

teta0=0;nh=1;ro=nh*sin(teta0*pi/180);
delta0=0;

nn=[4,4];% ordres de fourier 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= 1 ;   
textures{2}= inf  ; 
textures{3}={ inf   ,  [0,0,3,3,  1]     };
textures{4}={ inf   ,  [0,0,8,8,  1]     };

% initialisation
parm=res0;  % parametres par defaut
parm.sym.x=0;parm.sym.y=0;parm.sym.pol=1;% utilisation de 2 symetries
parm.res1.trace=1;%  trace des textures     (defaut 0)
parm.res1.champ=1; % si on veut un calcul soigne du champ

aa=res1(LD,D,textures,nn,ro,delta0,parm);


% recherche de la resonance en h comme pole de h-->R
er=inf;iter=0;z0=20;z=[];zz=[];
while(er>eps)&(iter<30);
profil={[0,1,z0,0] ,[1,3,4,2]  };
ef=res2(aa,profil,parm);
zz0=1/ef.TEinc_top_reflected.amplitude_TE{0};
iter=iter+1;z=[z;z0];zz=[zz;zz0];[z0,z,zz,er]=retcadillac(z,zz);
end;

h=real(z0);


rr=[];
hh=h+linspace(-.1,.1,30);
for hhh=hh;
profil={[0,1,hhh,0] ,[1,3,4,2]  };
ef=res2(aa,profil);
rr=[rr,ef.TEinc_top_reflected.amplitude_TE{0}];
end;
figure;plot(hh,real(rr),hh,imag(rr));xlabel('hauteur');legend('real R','imag R');title(['resonance pour h=',num2str(h)]);pause(eps);

% coupe y=0  
x=linspace(-1.5*D(1),1.5*D(1),150);y=0;
einc=[0,1];   %  composantes du champ e incident  TE

parm.res3.npts=[10,5,40,10];% nombre de points par tranche
profil={[5,1,h,5] , [1,3,4,2]  };
[e,z,o]=res3(x,y,aa,profil,einc,parm);
figure;

subplot(2,2,1);retcolor(x,z,real(o(:,:,1)));xlabel('Z');ylabel('X');axis equal;title(['objet coupe Y=',num2str(y)]);
subplot(2,2,2);retcolor(x,z,abs(e(:,:,1,2)).^2);xlabel('Z');ylabel('X');title('EY');axis equal;
subplot(2,2,3);retcolor(x,z,abs(e(:,:,1,4)).^2);xlabel('Z');ylabel('X');title('HX');axis equal;
subplot(2,2,4);retcolor(x,z,abs(e(:,:,1,6)).^2);xlabel('Z');ylabel('X');title('HZ');axis equal;


retio;