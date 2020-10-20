%   2  D    exemple1_2D
%  TEST DE CONVERGENCE POUR L'EFFACITEE D'UN RESEAU SIMPLE 


clear;
LD=6;% longueur d'onde
D=[10,10];% pas du reseau

teta0=0;nh=1;ro=nh*sin(teta0*pi/180);
delta0=0;
parm=res0;
parm.sym.x=0;parm.sym.y=0;parm.sym.pol=1;% utilisation de 2 symetries
figure;
t=[];ntf=[];
for n=0:10;%boucle sur les ordres de fourier
nn=[n,n]; 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= 1;   
textures{2}= 1.5 ; 
textures{3}={ 1.5   ,  [0,0,5,8,  1,  1]     };
% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);
% definition du profil et calcul de la diffraction
profil={[0,20,0], [1,3,2]  };

ef=res2(aa,profil);

% transmission dans l'ordre 0
t=[t,ef.TEinc_top_transmitted.efficiency_TE{0}];
ntf=[ntf,n];
plot(ntf,t,'.');xlabel('nombre de termes de fourier');ylabel('transmission dans l''ordre 0');title('TEST DE CONVERGENCE');pause(eps);
end

retio;



