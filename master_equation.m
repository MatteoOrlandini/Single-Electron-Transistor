%%% metodo master equations %%%

%V: tensione di drain
%Vg: tensione di gate
%N: elettroni presenti nel dot

function y = master_equation(V,Vg,N);

clear ma;
for i = -N : 1 : N
    %ma(N+i+1,N+i+2)  = - (gamma_L->dot + gamma_R->dot + gamma_dot->L + gamma_dot->R)
    ma(N+i+1,N+i+2) = -(f_tunnel(1,i,V,Vg) + f_tunnel(2,i,V,Vg) + f_tunnel(3,i,V,Vg) + f_tunnel(4,i,V,Vg));
    %ma(N+i+1,N+i+2) = gamma_L->dot + gamma_R->dot
    ma(N+i+1,N+i+2-1) = f_tunnel(1,i-1,V,Vg) + f_tunnel(4,i-1,V,Vg); 
    %ma(N+i+1,N+i+2+1) = gamma_dot->L è gamma_dot->R
    ma(N+i+1,N+i+2+1) = f_tunnel(2,i+1,V,Vg) + f_tunnel(3,i+1,V,Vg);
end;
mma = [1,zeros(1,2*N+2); ma; zeros(1,2*N+2),1];
%o è una matrice diagonale contenente gli autovalori, p è una matrice le
%cui colonne sono i corrispondenti autovettori (mma*p = p*o)
[p,o] = eig(mma);
%s è il minimo degli autovalori (contenuti nella diagonale della
%matrice o) e posiz è l'indice dell'autovalore nella diagonale
[s,posiz] = min(abs(diag(o)));
%p_n è la probabilità che n cariche sono presenti nel dot
%perchè si calcola così?????????????
p_n(:,1) = abs(p(:,posiz))/sum(abs(p(:,posiz)));
for i = 1:length(p_n) 
    %prodotto della probabilità p_n per la differenza
    % tra il tunneling rate dot->R e R->dot
    contr(i) = p_n(i)*(f_tunnel(1,i-N-2,V,Vg)-f_tunnel(2,i-N-2,V,Vg));
end;
%la somma di contr è la corrente
y = sum(contr);

