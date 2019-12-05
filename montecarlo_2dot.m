clear all; close all;

y = [];
n = 10;
N = -n:n;

for i = -n : n
    y = [y; f_tunnel_2dot(1,i,5,1,0)];
end
%%% plot dei 4 tipi di gamma(n) per n che va da -10 a 10
figure ('Name','plot dei 4 tipi di gamma(n) per n che va da -10 a 10','NumberTitle','off');
plot(N, y(:,1), N, y(:,2), N, y(:,3), N, y(:,4), N, y(:,5), N, y(:,6));
xlabel ('Numero di cariche nel dot (n)');
ylabel ('Rate di tunneling (\Gamma_n)');
legend ('\Gamma_{source->dot1}', '\Gamma_{dot1->source}', '\Gamma_{dot1->dot2}', '\Gamma_{dot2->dot1}', '\Gamma_{dot2->drain}', '\Gamma_{drain->dot2}')