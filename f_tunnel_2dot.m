% calcolo probabilità di tunneling per transistor con 2 dot

function y = f_tunnel_2dot(j,n1,n2,Vd,Vg)

e = 1;
Cg = 0.8;
Cd = 0.8;
Cs = 0.8;
Cdot = 0.8;
C = Cg+Cd+Cs+Cdot;

V1 = ((-n*e+Cg*Vg)*(C+Cdot) + Cdot*Cd*Vd)/(C^2-Cdot^2);
V2 = (-n*e+Cdot*V1+Cd*Vd+Cg*Vg)/C;

%notazione: dE_inizio_fine
dE_s_dot1 = -e*V1+(e^2)/2*C;
dE_dot1_s = e*V1+(e^2)/2*C;
dE_dot1_dot2 = e*V1-e*V2+2*(e^2)/2*C;
dE_dot2_dot1 = -e*V1+e*V2+2*(e^2)/2*C;
dE_dot2_d = -e*Vd+e*V2+(e^2)/2*C;
dE_d_dot2 = e*Vd-e*V2+(e^2)/2*C;

dE = [dE_s_dot1, dE_dot1_s, dE_dot1_dot2, dE_dot2_dot1, dE_dot2_d, dE_d_dot2];

%provvisorio 
y = [0.001*dE(1)/(exp(5*dE(1)-1)),0.001*dE(2)/(exp(5*dE(2)-1)),0.001*dE(3)/(exp(5*dE(3)-1)),0.001*dE(4)/(exp(5*dE(4)-1)),0.001*dE(5)/(exp(5*dE(5)-1)),0.001*dE(6)/(exp(5*dE(6)-1))];

