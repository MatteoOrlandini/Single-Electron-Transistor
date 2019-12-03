function y=f_tunnel0(j,n,V,Vg)


%%% n è il numero di cariche (elettroni) nel DOT (da -10 a 10 negli esempi)
%%% l'uscita è data dai 4 tipi di ?(n), infatti j va da 1 a 4, rispettivamente 1) da DOt a Source, 2) da Source a Dot, 3) da dot a Drain, 4) da Drain a Dot 

cg=0.8;cd=0.8;cs=0.8;
c=cg+cd+cs;
dE1=-(-n+cg*Vg+cd*V)/c+5/c;
dE2=(-n+cg*Vg+cd*V)/c+5/c;
dE3=-V+(-n+cg*Vg+cd*V)/c+5/c;
dE4=+V-(-n+cg*Vg+cd*V)/c+5/c;

f=0.001*[dE1/(exp(5*dE1)-1),dE2/(exp(5*dE2)-1),dE3/(exp(5*dE3)-1),dE4/(exp(5*dE4)-1)];


y=f; 
