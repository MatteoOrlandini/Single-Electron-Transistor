%%% Metodo Montecarlo %%%

%N: numero di elettroni presenti nel dot
%V: tensione di drain
%Vg: tensione di gate

function [tempo,q,s,d] = arrivi(N,V,Vg);

tempo = zeros(N,1);
%q = elettroni del quantum dot
q = zeros(N,1);
%s = elettroni nel source
s = zeros(N,1);
%d = elettroni nel drain
d = zeros(N,1);
%N sono le iterazioni del metodo montecarlo (di solito 1000)
for i = 1:N
    for j = 1:4
        %r(j) è un numero casuale
        r(j) = rand(1,1);
        %t(j) è una matrice che contiene il tempo
        t(j) = -log(r(j))/f_tunnel(j,q(i),V,Vg);
    end;
    %k = tempo minimo, h = indice del tempo minimo da cui si capisce la
    %gamma minima
    [k,h] = min(t);
    tempo(i) = k;
    %se h = 1, è il caso in cui un elettrone va dal source al dot: il
    %source perde un elettrone, il dot acquista un elettrone, 
    %il drain rimane uguale
    if h == 1
        s(i+1) = s(i)-1;
        q(i+1) = q(i)+1;
        d(i+1) = d(i);
    %se h == length(t) (== 4), è il caso in cui un elettrone da dal drain
    %al dot: il drain perde un elettrone, il dot ne acquista uno e il
    %source rimane uguale
    else if h == length(t)
            d(i+1) = d(i)-1;
            q(i+1) = q(i)+1;
            s(i+1) = s(i);
        %se h = 2, è il caso in cui un elettrone va dal dot al source: il
        %source acquista un elettrone, il dot perde un elettrone, 
        %il drain rimane uguale
        else if h == 2
                s(i+1) = s(i)+1;
                q(i+1) = q(i)-1;
                d(i+1) = d(i);
            %il caso rimanente (h == 3) è quello in cui in cui un elettrone
            %va dal dot al drain: il source rimane uguale, il dot perde un 
            %elettrone, il drain acquista un elettrone
            else
                d(i+1) = d(i)+1;
                q(i+1) = q(i)-1;  
                s(i+1) = s(i);
            end;
        end;
    end;
end;




