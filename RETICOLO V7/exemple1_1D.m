%   1  D    exemple1_1D
%  TEST DE CONVERGENCE POUR L'EFFACITEE D'UN RESEAU SIMPLE 



clear;
LD=6;% longueur d'onde
D=10;% pas du reseau
pol=-1; % 1:TE   -1:TM
%pol=1; % 1:TE   -1:TM
parm.sym.x=0;% utilisation de la symetrie
figure;
t=[];ntf=[];
for n=5:2:50;%boucle sur les ordres de fourier
    
nh=1;teta0=0;beta0=nh*sin(teta0*pi/180);
parm=res0(pol);% initialisation des parametres par defaut

parm.sym.x=0;% utilisation de la symetrie

 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures{1}= nh;   
textures{2}= 1.5 ; 
textures{3}={[-4,4],[1.5,1]  };
% initialisation
aa=res1(LD,D,textures,n,beta0,parm);
% definition du profil et calcul de la diffraction
profil={[0,20,0], [1,3,2]  };

ef=res2(aa,profil);

% recherche d'un ordre


% efficacite de l'onde  diffractee en bas dans l'ordre 0 pour une onde incidente du haut
t=[t,ef.inc_bottom_transmitted.efficiency{0}];
ntf=[ntf,n];
plot(ntf,t,'.');xlabel('nombre de termes de fourier');ylabel('transmission dans l''ordre 0');grid;title('TEST DE CONVERGENCE');pause(eps);
end

retio;



