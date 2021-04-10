function tunneling_rate = f_tunnel0(n, Vd, Vg)
%f_tunnel0 computes the tunneling rate 
%   tunneling_rate = f_tunnel0(n, Vd, Vg) stores the tunneling_rate  in a 
%   matrix for n electrons, drain voltage Vd and gate voltage Vg.
%   The first column of tunneling_rate is the tunneling rate from dot to 
%   source, the second column is from source to dot, the third column is 
%   from dot to drain and the fourth is from drain to dot

cg = 0.8; % gate capacitance
cd = 0.8; % drain capacitance
cs = 0.8; % source capacitance
c = cg+cd+cs; % total capacitance
% computing differential energies
differentialEnergy1 = -(-n+cg*Vg+cd*Vd)/c+5/c;    % from dot to source
differentialEnergy2 = (-n+cg*Vg+cd*Vd)/c+5/c;     % from source to dot
differentialEnergy3 = -Vd+(-n+cg*Vg+cd*Vd)/c+5/c; % from dot to drain
differentialEnergy4 = +Vd-(-n+cg*Vg+cd*Vd)/c+5/c; % from drain to dot

tunneling_rate = 0.001*[differentialEnergy1/(exp(5*differentialEnergy1)-1)...
    differentialEnergy2/(exp(5*differentialEnergy2)-1)...
    differentialEnergy3/(exp(5*differentialEnergy3)-1)...
    differentialEnergy4/(exp(5*differentialEnergy4)-1)];