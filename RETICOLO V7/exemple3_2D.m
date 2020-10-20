%   2  D    exemple3_2D
%  TRACE DU CHAMP  2D


clear;
LD=6;% longueur d'onde
D=[10,10];% pas du reseau

teta0=0;nh=1;ro=nh*sin(teta0*pi/180);
delta0=0;
parm=res0;
parm.sym.x=0;parm.sym.y=0;parm.sym.pol=1;% utilisation de 2 symetries
parm.res1.champ=1; % si on veut un calcul soigne du champ

nn=[5,5];% ordres de fourier 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}=1;   
textures{2}= 1.5; 
textures{3}={ 1.5   ,  [0,0,5,8,  1,  1]     };
% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);

profil={[3,20,3] ,[1,3,2]  };
    

% coupe y=0  
x=linspace(0,D(1),150);y=0;
einc=[0,1];   %  composantes du champ e incident  TE
parm=res0; %parametres par defaut
parm.res3.npts=[5,40,5];% nombre de points par tranche

[e,z,o]=res3(x,y,aa,profil,einc,parm);

figure;

retsubplot(1,4,1);retcolor(x,z,real(o(:,:,1)));xlabel('Z');ylabel('X');axis equal;title(['objet coupe Y=',num2str(y)]);
retsubplot(1,4,2);retcolor(x,z,abs(e(:,:,1,2)).^2);xlabel('Z');ylabel('X');title('EY');axis equal;
retsubplot(1,4,3);retcolor(x,z,abs(e(:,:,1,4)).^2);xlabel('Z');ylabel('X');title('HX');axis equal;
retsubplot(1,4,4);retcolor(x,z,abs(e(:,:,1,6)).^2);xlabel('Z');ylabel('X');title('HZ');axis equal;

retio;