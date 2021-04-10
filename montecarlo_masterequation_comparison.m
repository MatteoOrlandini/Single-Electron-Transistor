clc
clear all
close all

%% Tunneling rate
n = 10; % number of charges stored into the dot
Vd = 1; % drain voltage
Vg = 0; % gate voltage
y = zeros(2*n+1,4);
for i = -n : n 
    y(i+n+1,:) = f_tunnel0(i,Vd,Vg);
end
% plot gamma(n) for n from -10 to 10
figure ('Name', "Tunneling rate gamma(n) for n from -"+n+" to "+n, 'NumberTitle', 'off');
plot(-n:n, y(:,1), -n:n, y(:,2), -n:n, y(:,3), -n:n, y(:,4));
xlabel ('Number of charges stored into the dot (n)', 'Interpreter', 'latex');
ylabel ('Tunneling rate $\Gamma_n$',  'Interpreter', 'latex');
title("V_d = " + Vd + " V, V_g = " + Vg + " V , " + -n + " < n < " + n);
legend ('\Gamma_{source \rightarrow dot}', '\Gamma_{dot \rightarrow source}', '\Gamma_{dot \rightarrow drain}', '\Gamma_{drain \rightarrow dot}', 'Location', 'best')

%% Master equation - Computing current density and density
% plot p(n)*[(Gamma_dot->R)(N) - (Gamma_R->dot)(N)] vs. Vd for different Vg
n = 10; % number of charges stored into the dot
Vg = 0.5; % gate voltage
current = zeros(50, 1); % calculate 50 values of current
figure('Name', 'Charge density for different drain voltages', 'NumberTitle','off');
for i = 1:50
    Vd = exp(0.04*i)-1; % drain voltage
    [charge_density, current(i)] = master_equation0(Vd, Vg, n); 
    plot(-n-1: n+1, charge_density)
    hold on
end
xlabel('Number of electrons stored into the dot [n]');
xlim([-n-1, n+1]);
ylabel ('$p(n) \cdot [\Gamma_{dot \rightarrow R}(N) - \Gamma_{R \rightarrow dot}(N)]$', 'Interpreter', 'latex');
title("V_g = " + Vg + " V, n = " + n);

figure('Name', 'Current vs. drain voltage', 'NumberTitle','off');
plot(exp(0.04*[1:50])-1, current)
xlabel ('Drain voltage ($V_d$)', 'Interpreter', 'latex');
ylabel ('Current');
title("V_g = " + Vg + " V, n = " + n);

% plot  p(n)*[(Gamma_dot->R)(N) - (Gamma_R->dot)(N)] vs. Vg for different Vd
figure('Name', 'Charge density for different gate voltages', 'NumberTitle','off');
Vd = 0.1; % drain voltage
for i=1:50 
    Vg = exp(0.04*i)-1; % gate voltage
    [charge_density, current(i)] = master_equation0(Vd, Vg, n);
    plot(-n-1: n+1, charge_density)
    hold on
end
xlabel('Number of electrons stored into the dot [n]');
xlim([-n-1, n+1]);
ylabel ('$p(n) \cdot [\Gamma_{dot \rightarrow R}(N) - \Gamma_{R \rightarrow dot}(N)]$', 'Interpreter', 'latex');
title("V_d = " + Vd + " V, n = " + n);

figure('Name', 'Current vs. gate voltage', 'NumberTitle','off');
plot(exp(0.04*[1:50])-1, current)
xlabel ('Gate voltage ($V_g$)', 'Interpreter', 'latex');
ylabel ('Current');
title("V_d = " + Vd + " V, n = " + n);
%% Comparison between master equation and Monte Carlo, gate voltage fixed
clear all;
n = 10;   % number of electrons stored into the dot
Vg = 0.5; % gate voltage
current = zeros(50, 1); % calculate 50 values of current
for i = 1:50 
    Vd = exp(0.04*i)-1; % drain voltage
    current(i) = master_equation(Vd, Vg, n); % compute master equation
end
N = 10000;  % number of Monte Carlo iterations
Vg = 0.5;   % gate voltage
for i = 1:50 
    Vd = exp(0.04*i)-1; % drain voltage
    [time,q,s,d] = monte_carlo(N, Vd, Vg); % compute Monte Carlo method
    drain_current(:,i) = d(:,1);     % drain current
    source_current(:,i) = s(:,1);    % source current
    t(:,i) = time(:,1);              % time
    dot_charge(:,i) = q(:,1);        % dot charge
end
% plot Drain current vs. Vd, Vg is fixed, Vg=0.5, comparison with master equation result
figure ('Name', 'Drain current vs. Vd', 'NumberTitle','off');
plot(exp(0.04*[1:50])-1, (drain_current(length(drain_current),:))./sum(t(:,:)),...
    exp(0.04*[1:50])-1, current);
xlabel ('Drain voltage ($V_d$)', 'Interpreter', 'latex');
ylabel ('Drain current ($I_d$)', 'Interpreter', 'latex');
title("Drain current vs. V_d, V_g = " + Vg + " V, comparison with master equation result");
legend ('Monte Carlo', 'Master Equations', 'Location', 'best');

for i = 1:N 
    partial_current(i) = drain_current(i,10)/sum(t(1:i,10));
end
% plot Drain current vs. Time, Vd and Vg are fixed, Vd=exp(0.04*10), Vg=0.5
figure ('Name', 'Drain current vs. Time', 'NumberTitle', 'off');
plot(partial_current);
xlabel ('Time');
xlim([0 length(partial_current)])
ylabel ('Drain current ($I_d$)', 'Interpreter', 'latex');
title("Drain current vs. Time, V_d = " + exp(0.04*10) + " V, V_g = " + Vg + " V");

%%% plot number of charges arrived at the drain vs. Time, Vd and Vg are fixed, Vd=exp(0.04*10), Vg=0.5
figure ('Name', 'Number of charges arrived at the drain vs. Time', 'NumberTitle','off');
plot(drain_current(:,10));
xlabel ('Time');
xlim([0 length(drain_current)])
ylabel ('Charges arrived at the drain');
title("Number of charges arrived at the drain vs. Time, V_d = " + exp(0.04*10) + " V, V_g = " + Vg + " V");

%%% plot number of charges that left the drain vs. Time, Vd and Vg are fixed, Vd=exp(0.04*10), Vg=0.5
figure ('Name', 'Number of charges that left the drain vs. Time', 'NumberTitle','off');
plot(source_current(:,10))
xlabel ('Time');
xlim([0 length(source_current)])
ylabel ('Charges that left the drain');
title("Number of charges that left the drain vs. Time, V_d = " + exp(0.04*10) + " V, V_g = " + Vg + " V");

%%% plot algebraic sum of the number of charges referred to in the two previous points
figure ('Name', 'Algebraic sum of the number of charges arrived at the drain and that left the drain', 'NumberTitle','off');
plot(source_current(:,10) + drain_current(:,10))
xlabel ('Time');
xlim([0 length(source_current)])
ylabel ('Charges that left the source');
title("Algebraic sum of the number of charges arrived at the drain and that left the drain, V_d = " + exp(0.04*10) + " V, V_g = " + Vg + " V");

% temporal average of the algebraic sum referred to in the previous point (sum of the number of charges arrived at the drain and that left the drain)
disp ("Temporal average of the number of charges that arrived and left the drain: " + ...
    mean(drain_current(:,10)+source_current(:,10)));
% temporal average of the charge stored into the dot, Vd and Vg are fixed, Vd = exp(0.04 * 10) and Vg = 0.5 V
disp ("Temporal average of the number of charges stored into the dot: " + ...
    mean(dot_charge(:,10)));
disp ("Note the equality of the two averages referred to in the two previous points (sign aside, being a balance of offices)")

%% Comparison between master equation and Monte Carlo, fixed drain voltage 
clear all; 
n = 10; % number of electrons
Vd = 0.4;    % drain voltage
for i = 1:50 
    Vg = exp(0.04*i)-1; % gate voltage
    current(i) = master_equation(Vd, Vg, n); 
end
N = 10000;  % number of Monte Carlo iterations
for i = 1:50 
    Vg = exp(0.04*i)-1;
    [time, q, s, d] = monte_carlo(N, Vd, Vg); 
    drain_current(:,i) = d(:);     % drain current
    source_current(:,i) = s(:,1);  % source current
    t(:,i) = time(:,1);            % time
    dot_charge(:,i) = q(:,1);      % total charge
end

% plot drain current vs. Vg, Vd is fixed, Vd=0.4, comparison with master equation result
figure ('Name', 'Drain current vs. Vg, comparison with master equation result', 'NumberTitle','off');
plot(exp(0.04*[1:50])-1, (drain_current(length(drain_current),:))./sum(t(:,:)),...
    exp(0.04*[1:50])-1, current);
xlabel ('Gate voltage ($V_g$)', 'Interpreter', 'latex');
ylabel ('Drain current ($I_d$)', 'Interpreter', 'latex');
title ("Drain current vs. V_g, V_d = " + Vd + " V, comparison with master equation result");
legend ('Monte Carlo', 'Master Equations');

for i = 1:N 
    partial_current(i) = drain_current(i,10)/sum(t(1:i,10));
end

% plot Drain current vs. Time, Vg and Vd are fixed, Vg=exp(0.04*10) and Vd=0.4
figure ('Name', 'Drain current vs. Time','NumberTitle','off');
plot(partial_current)  
xlabel ('Time');
ylabel ('Drain current ($I_d$)', 'Interpreter', 'latex');
title ("Drain current vs. Time, V_g = " + exp(0.04*10) + " V, V_d = " + Vd);

%% Comparison between master equation and Monte Carlo varying drain and gate voltage
clear all
n = 10; % number of electrons stored in dot
for i = 1:50 
    Vg = exp(0.04*i)-1; % gate voltage
    for j = 1:20 
        Vd = exp(0.04*j)-1; % drain voltage
        current(i,j) = master_equation(Vd, Vg, n); % compute master equation method
    end
end
N = 10000; % number of Monte Carlo iterations
for i = 1:50 
    Vg = exp(0.04*i)-1; % gate voltage
    for j = 1:20 
        Vd = exp(0.04*j)-1; % drain voltage
        [time, q, s, d] = monte_carlo(N, Vd, Vg); % compute Monte Carlo method
        drain_current(i,j) = d(length(d))/sum(time);    % drain current
        source_current(i,j) = s(length(s))/sum(time);   % source current
        dot_charge(i,j) = mean(q);  % dot charge
    end
end
% plot comparison Monte Carlo and Master equations results
figure('Name', 'Comparison Monte Carlo and Master equations results', 'NumberTitle','off');
% master equation results
surf(exp(0.04*[1:20])-1, exp(0.04*[1:50])-1, current, 'FaceColor', 'red'); 
xlabel('Drain voltage $V_d$', 'Interpreter', 'latex');
ylabel('Gate voltage $V_g$', 'Interpreter', 'latex');
hold on; 
% montecarlo results
surf(exp(0.04*[1:20])-1, exp(0.04*[1:50])-1, drain_current, 'FaceColor', 'blue'); 
xlabel('Drain voltage $V_d$', 'Interpreter', 'latex');
ylabel('Gate voltage $V_g$', 'Interpreter', 'latex');
zlabel ('Drain current');
title ('Comparison Monte Carlo and Master equations results, varying V_d and V_g');
legend('Master Equation results', 'Monte Carlo results', 'Location', 'eastoutside');

% plot charges stored into the dot, varying Vd and Vg
figure('Name', 'Charges stored into the dot, varying Vd and Vg', 'NumberTitle',' off');
surf(exp(0.04*[1:20])-1, exp(0.04*[1:50])-1, dot_charge);
xlabel('Drain voltage $V_d$', 'Interpreter', 'latex');
ylabel('Gate voltage $V_g$', 'Interpreter', 'latex');
zlabel('Charges stored into the dot');
title('Charges stored into the dot, varying V_d and V_g');
c = colorbar('Location', 'eastoutside');
c.Label.String = 'Charges stored into the dot';