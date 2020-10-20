%   2  D    exemple1_conique
%  INCIDENCE CONIQUE 

clear;
LD=6;% longueur d'onde
D=10;% pas du reseau

teta0=30;nh=1;ro=nh*sin(teta0*pi/180);
delta0=20;
parm=res0;

parm.res1.champ=1; % si on veut un calcul soigne du champ

nn=5;% ordres de fourier 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= 1;   
textures{2}= 1.5; 
textures{3}={ [-2.5,2.5] , [1.5,1]};
% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);
% definition du profil et calcul de la diffraction 
profil={[0,20,0] ,[1,3,2]  };

[ef,tab]=res2(aa,profil);
    

%efficacite transmise et reflechie 
T=ef.TEinc_top_transmitted.efficiency_TE{0};
R=ef.TEinc_top_reflected.efficiency_TE{0};
disp(rettexte(R,T));

% coupe y=0  
x=linspace(0,D(1),150);

parm=res0;parm.res2.result=0;            %parametres par defaut
parm.res3.npts=[5,40,5];% nombre de points par tranche
profil={[5,20,5] ,[1,3,2]  };

einc=[0,1]; [eTE,z,o]=res3(x,aa,profil,einc,parm);  %   incident  TE
einc=[1,0]; [eTM,z,o]=res3(x,aa,profil,einc,parm);  %   incident  TE

figure;
subplot(3,3,2);retcolor(x,z,real(o));title('objet coupe ');

for psi=0:20:180;e=eTE*cos(psi*pi/180)+eTM*sin(psi*pi/180);
subplot(3,3,3);plot([0,cos(psi*pi/180)],[0,sin(psi*pi/180)]);axis([-1,1,-1,1]);title('polarisation');    
subplot(3,3,4);retcolor(x,z,abs(e(:,:,1)).^2);title('EX');    
subplot(3,3,5);retcolor(x,z,abs(e(:,:,2)).^2);title('EY');
subplot(3,3,6);retcolor(x,z,abs(e(:,:,3)).^2);title('EZ');
subplot(3,3,7);retcolor(x,z,abs(e(:,:,4)).^2);;title('HX');
subplot(3,3,8);retcolor(x,z,abs(e(:,:,5)).^2);title('HY');
subplot(3,3,9);retcolor(x,z,abs(e(:,:,6)).^2);;title('HZ');
drawnow;
end;
retio;

