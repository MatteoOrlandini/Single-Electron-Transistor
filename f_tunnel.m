function tunneling_rate = f_tunnel(j,n,V,Vg)
%f_tunnel computes the tunneling rate
%   tunneling_rate = f_tunnel(j,n,V,Vg) returns a scalar tunneling_rate  
%   for n electrons, drain voltage Vd and gate voltage Vg. 
%   If j = 1, f_tunnel calculates the tunneling rate from dot to source
%   If j = 2, f_tunnel calculates the tunneling rate from source to dot
%   If j = 3, f_tunnel calculates the tunneling rate from dot to drain
%   If j = 4, f_tunnel calculates the tunneling rate from drain to dot

cg = 0.8; % gate capacitance
cd = 0.8; % drain capacitance
cs = 0.8; % source capacitance
c = cg+cd+cs; % total capacitance
c = cg+cd+cs;
%delta energy source->dot
dE1 = -(-n+cg*Vg+cd*V)/c+5/c;
%delta energy dot->source
dE2 = (-n+cg*Vg+cd*V)/c+5/c;
%delta energy dot->drain
dE3 = -V+(-n+cg*Vg+cd*V)/c+5/c;
%delta energy drain->dot
dE4 = +V-(-n+cg*Vg+cd*V)/c+5/c;

% tunneling rate array
f = 0.001*[dE1/(exp(5*dE1)-1), dE2/(exp(5*dE2)-1), dE3/(exp(5*dE3)-1), dE4/(exp(5*dE4)-1)];

%j = 1: gamma_source->dot
%j = 2: gamma_dot->source
%j = 3: gamma_dot->drain
%j = 4: gamma_drain->dot

tunneling_rate = f(j);
