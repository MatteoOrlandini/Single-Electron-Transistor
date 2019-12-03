%%% calcolo della probabilità di tunneling %%%

%j: indice della j-esima gamma
%j = 1: gamma_source->dot
%j = 2: gamma_dot->source
%j = 3: gamma_dot->drain
%j = 4: gamma_drain->dot
%n: numero di elettroni presenti nel dot
%V: tensione di drain
%Vg: tensione di gate

function y = f_tunnel(j,n,V,Vg)

% n è il numero di cariche (elettroni) nel DOT (da -10 a 10 negli esempi)
% l'uscita è data dai 4 tipi di ?(n), infatti j va da 1 a 4, rispettivamente 
% 1) da Dot a Source, 2) da Source a Dot, 3) da Dot a Drain, 4) da Drain a Dot 

cg = 0.8;
cd = 0.8;
cs = 0.8;
c = cg+cd+cs;
%delta E source->dot
dE1 = -(-n+cg*Vg+cd*V)/c+5/c;
%delta E dot->source
dE2 = (-n+cg*Vg+cd*V)/c+5/c;
%delta E dot->drain
dE3 = -V+(-n+cg*Vg+cd*V)/c+5/c;
%delta E drain->dot
dE4 = +V-(-n+cg*Vg+cd*V)/c+5/c;

%f = rate di tunneling (gamma_i)
f = 0.001*[dE1/(exp(5*dE1)-1), dE2/(exp(5*dE2)-1), dE3/(exp(5*dE3)-1), dE4/(exp(5*dE4)-1)];
%f(j) contiene il rate di tunneling di una delle quattro possibilità di
%tunneling
y=f(j);

