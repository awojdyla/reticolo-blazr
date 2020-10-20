%   2  D    exemple10_2D
%  RESEAU ECHELETTE 1D cas particulier de l'incidence conique avec delta=0    
%  (comparer avec exemple10_1D)

clear;

te=[];tm=[];teta=[];
premier=1;

for teta0=linspace(-30,89,91);
LD=13;% longueur d'onde
D=30;% pas du reseau
if teta0>0;delta0=0;else delta0=180;end;nh=1;ro=nh*sin(abs(teta0)*pi/180);% car dans le repere des ondes planes teta>0
for pol=[1,-1];
parm=res0;
parm.sym.y=0;%  symetries par rapport a y=0 
parm.sym.pol=pol;  % 1 TE  -1 TM


nn=5;% ordres de fourier 

% description des textures
for ii=1:12;textures{ii}={ [-ii*D/13,0] , [1,1.5]};end;
textures{13}= 1;   
textures{14}= 1.5 ; 

% initialisation
aa=res1(LD,D,textures,nn,ro,delta0,parm);
% definition du profil et calcul de la diffraction 
profil={[0,ones(1,12)*(D/13),0],  [13,1:12,14]  };

%au premier passage on trace une coupe du reseau (pour verifier)
if premier==1;
figure;    
x=linspace(-D,D,100);    
parm.res3.cale=[] ;%signifie que l'on ne calcule pas le champ
profil1={[5,ones(1,12)*(D/13),5],  [13,1:12,14] };
[tab1,z,o]=res3(x,aa,profil1,[1,1],parm); % profil1 
retcolor(x,z,real(o));xlabel('X');ylabel('Z');title('coupe profil');axis equal;pause(eps)
figure;
premier=0;    
end;    


[ef,tab]=res2(aa,profil);
    


if pol==1;
te=[te,ef.TEinc_top_transmitted.efficiency_TE{1}];% TE incident en haut  transmission TE dans l'ordre 0
teta=[teta,teta0];
else;
tm=[tm,ef.TMinc_top_transmitted.efficiency_TM{1}];% TE incident en haut  transmission TE dans l'ordre 0
plot(teta,te,teta,tm,'--');xlabel('teta');title('RESEAU 1D');legend('TE','TM');ylabel('efficacitee transmise dans l''ordre 1');pause(eps);
end;

end;

end;
retio
