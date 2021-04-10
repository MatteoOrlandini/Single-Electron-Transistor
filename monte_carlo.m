function [time, q, s, d] = monte_carlo(n, Vd, Vg)
%monte_carlo computes the Monte Carlo method
%   [time, q, s, d] = monte_carlo(n, Vd, Vg) returns "silence time" array, 
%   q electrons stored into the quantum dot, s electrons in source, d 
%   electrons in drain for n iterations, drain voltage Vd and gate voltage Vg.

% time array
time = zeros(n,1);
% q: electrons quantum dot
q = zeros(n,1);
% s: electrons in source
s = zeros(n,1);
% d: electrons in drain
d = zeros(n,1);

for i = 1:n   % n: Monte Carlo iterations
    for j = 1:4 
        r(j) = rand(1,1);   %r(j) random number
        t(j) = -log(r(j))/f_tunnel(j, q(i), Vd, Vg); %t(j) "silence time" matrix
    end
    [k,h] = min(t); % k minimum time, h index of the minimum "silencetime" (required for minimum gamma)
    time(i) = k;
    
    if (h == 1)   % electron goes from source to dot
        s(i+1) = s(i)-1;    % source losts an electron
        q(i+1) = q(i)+1;    % dot gets an electron
        d(i+1) = d(i);      % drain electrons are same as the previous step 
        
    elseif (h == length(t))  % electron goes from drain to dot (length(t) == 4)
            d(i+1) = d(i)-1;    % drain losts an electron
            q(i+1) = q(i)+1;    % dot gets an electron
            s(i+1) = s(i);      % source electrons are same as the previous step 

    elseif (h == 2)  % electron goes from dot to source
            s(i+1) = s(i)+1;    % source gets an electron
            q(i+1) = q(i)-1;    % dot losts an electron
            d(i+1) = d(i);      % drain electrons are same as the previous step 
            
    elseif (h == 3)     % electron goes from dot to drain
        d(i+1) = d(i)+1;    % drain gets an electron
        q(i+1) = q(i)-1;    % dot losts an electron
        s(i+1) = s(i);      % source electrons are same as the previous step 
    end
end
end




