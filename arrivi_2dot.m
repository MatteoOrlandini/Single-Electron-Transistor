%%% Metodo Montecarlo %%%

%N: numero iterazioni
%Vd: tensione di drain
%Vg: tensione di gate

function [time,q1,q2,s,d] = arrivi_2dot(N,Vd,Vg)

time = zeros(N,1);
%q1 = elettroni del primo quantum dot
q1 = zeros(N,1);
%q2 = elettroni del secondo quantum dot
q2 = zeros(N,1);
%s = elettroni nel source
s = zeros(N,1);
%d = elettroni nel drain
d = zeros(N,1);
%N sono le iterazioni del metodo montecarlo (di solito 1000)
for i = 1:N
    for j = 1:6
        %r(j) è un numero casuale
        r(j) = rand(1,1);
        %t(j) è una matrice che contiene il tempo
        t(j) = -log(r(j))/f_tunnel_2dot(j,q1(i),q2(i),Vd,Vg);
    end;
    %k = tempo minimo, h = indice del tempo minimo da cui si capisce la
    %gamma minima
    [k,h] = min(t);
    time(i) = k;
    %se h = 1, è il caso in cui un elettrone va dal source al dot1: il
    %source perde un elettrone, il dot1 acquista un elettrone,
    %il drain e il dot2 rimangono uguali
    if h == 1
        s(i+1) = s(i)-1;
        q1(i+1) = q1(i)+1;
        q2(i+1) = q2(i);
        d(i+1) = d(i);
        %se h = 2, è il caso in cui un elettrone va dal dot1
        %al source
    else if h == 2
            s(i+1) = s(i)+1;
            q1(i+1) = q1(i)-1;
            q2(i+1) = q2(i);
            d(i+1) = d(i);
        else if h == 3
                s(i+1) = s(i);
                q1(i+1) = q1(i)-1;
                q2(i+1) = q2(i)+1;
                d(i+1) = d(i);
            else if h == 4
                    s(i+1) = s(i);
                    q1(i+1) = q1(i)+1;
                    q2(i+1) = q2(i)-1;
                    d(i+1) = d(i);
                else if h == 5
                        s(i+1) = s(i);
                        q1(i+1) = q1(i);
                        q2(i+1) = q2(i)-1;
                        d(i+1) = d(i)+1;
                    else if h == 6
                            s(i+1) = s(i);
                            q1(i+1) = q1(i);
                            q2(i+1) = q2(i)+1;
                            d(i+1) = d(i)-1;
                        end
                    end
                end
            end
        end
    end
end    
    
    
    
