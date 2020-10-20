%   2  D   exemple2_2D
%  VARIATION DE LA TRANSMISSION  FONCTION DE L'EPAISSEUR

clear;
LD=6;% longueur d'onde
D=[10,10];% pas du reseau

teta0=0;nh=1;ro=nh*sin(teta0*pi/180);
delta0=0;
parm=res0;
parm.sym.x=0;parm.sym.y=0;parm.sym.pol=1;% utilisation de 2 symetries
%parm.res1.trace=1; % trace des textures


nn=[5,5];% ordres de fourier 

% description des textures
textures{1}= 1;   
textures{2}= 1.5  ; 
textures{3}={ 1.5   ,  [0,0,5,8,  1,  1]     };
% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);
% definition du profil et calcul de la diffraction (on ne fait ce calcul qu'une fois)
figure;
t=[];epaisseur=[];
for h=linspace(0,40,50);

profil={[0,h,0] ,[1,3,2]  };
ef=res2(aa,profil);
 % TE diffractee en bas dans l'ordre 0 
t=[t,ef.TEinc_top_transmitted.efficiency{0}];
epaisseur=[epaisseur,h];
plot(epaisseur,abs(t).^2);xlabel('epaisseur');ylabel('transmission dans l''ordre 0');title('VARIATION DE LA TRANSMISSION DANS L''ORDRE 0 FONCTION DE L''EPAISSEUR');pause(eps);
end

retio;



