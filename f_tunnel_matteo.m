% calcolo probabilità di tunneling per transistor con 2 dot

function y = f_tunnel_matteo(n1,n2,Vd,Vg)
%definizione delle costanti 
e = 1.602*10^(-19);    % electronic charge (C)
atto = 10^(-18);    
Cg = 1 * atto;          %gate1 capacitor (F)
Cg1 = Cg/2 ;          %gate2 capacitor (F)
Cg2 = Cg/2;
Cd = 1 * atto;        %drain capacitor (F)
Cs = 1 * atto;        %source capacitor (F)
Cdot = 1 * atto;        %dot capacitor (dot1-dot2) (F)
Ctot = Cg+Cd+Cs+Cdot; %total capacitance (F)
temp = 10; % temperature
kb = 1.381*10^(-23);   % Boltzman constant (J/K)
mega = 10^6;     
r1 = 300*mega; % tunnel resistance R1 (Ohm)
r2 = 200*mega; % tunnel resistance R2 (Ohm)
r3 = 100*mega; % tunnel resistance R3 (Ohm)
Vg2=Vg/2;
Vg1=Vg/2;

%formule per due gate e due capacità di gate
%syms n1 e Vg1 Cg1 Vdot2 Cdot Cs n2 Vdot1 Vg2 Cg2 Vd Cd
%eqn1 = Vdot1 == (-n1*e + Vg1*Cg1 + Vdot2*Cdot)/(Cs + Cg1 + Cdot); %tensione dot 1
%eqn2 = Vdot2 ==(-n2*e + Vdot1*Cdot + Vg2*Cg2 + Vd*Cd)/(Cdot + Cg2 + Cd) ; %tensione dot 2
%y = solve(eqn1,eqn2,Vdot1,Vdot2); 

%formule per un gate e una capacità di gate
syms n1 e Vg1 Vdot2 Cdot Cs n2 Vdot1 Vg Vd Cd
eqn1 = Vdot1 == (-n1*e + Vg*Cg + Vdot2*Cdot)/(Cs + Cg + Cdot); %tensione dot 1
eqn2 = Vdot2 ==(-n2*e + Vdot1*Cdot + Vg*Cg + Vd*Cd)/(Cdot + Cg + Cd) ; %tensione dot 2
y = solve(eqn1,eqn2,Vdot1,Vdot2); 


Vdot1 = (Cd*Cdot*Vd + Cd*Cg1*Vg1 + Cdot*Cg1*Vg1 + Cdot*Cg2*Vg2 + Cg1*Cg2*Vg1 - Cd*e*n1 - Cdot*e*n1 - Cdot*e*n2 - Cg2*e*n1)/(Cd*Cdot + Cd*Cg1 + Cdot*Cg1 + Cdot*Cg2 + Cg1*Cg2 + Cd*Cs + Cdot*Cs + Cg2*Cs);
Vdot2 = (Cd*Cdot*Vd + Cd*Cg1*Vd + Cdot*Cg1*Vg1 + Cdot*Cg2*Vg2 + Cg1*Cg2*Vg2 + Cd*Cs*Vd + Cg2*Cs*Vg2 - Cdot*e*n1 - Cdot*e*n2 - Cg1*e*n2 - Cs*e*n2)/(Cd*Cdot + Cd*Cg1 + Cdot*Cg1 + Cdot*Cg2 + Cg1*Cg2 + Cd*Cs + Cdot*Cs + Cg2*Cs);


%notazione: dE_inizio_fine
dE_s_dot1 = -e*Vdot1 + (e^2)/2*Ctot;
dE_dot1_s = e*Vdot1 + (e^2)/2*Ctot;
dE_dot1_dot2 = e*Vdot1 - e*Vdot2+2*(e^2)/2*Ctot;
dE_dot2_dot1 = -e*Vdot1 + e*Vdot2+2*(e^2)/2*Ctot;
dE_dot2_d = -e*Vd + e*Vdot2+(e^2)/2*Ctot;
dE_d_dot2 = e*Vd - e*Vdot2+(e^2)/2*Ctot;

%dE = [dE_s_dot1, dE_dot1_s, dE_dot1_dot2, dE_dot2_dot1, dE_dot2_d, dE_d_dot2];

%provvisorio 
%y = [0.001*dE(1)/(exp(5*dE(1)-1)),0.001*dE(2)/(exp(5*dE(2)-1)),0.001*dE(3)/(exp(5*dE(3)-1)),0.001*dE(4)/(exp(5*dE(4)-1)),0.001*dE(5)/(exp(5*dE(5)-1)),0.001*dE(6)/(exp(5*dE(6)-1))];
%
gamma_s_dot1 = 1/(r1*e*e)*(-dE_s_dot1)/(1-exp(dE_s_dot1/(kb*temp)));
gamma_dot1_s = 1/(r1*e*e)*(-dE_dot1_s)/(1-exp(dE_dot1_s/(kb*temp)));
gamma_dot1_dot2 = 1/(r2*e*e)*(-dE_dot1_dot2)/(1-exp(dE_dot1_dot2/(kb*temp)));
gamma_dot2_dot1 = 1/(r2*e*e)*(-dE_dot2_dot1)/(1-exp(dE_dot2_dot1/(kb*temp)));
gamma_dot2_d = 1/(r3*e*e)*(-dE_dot2_d)/(1-exp(dE_dot2_d/(kb*temp)));
gamma_d_dot2 = 1/(r3*e*e)*(-dE_d_dot2)/(1-exp(dE_d_dot2/(kb*temp)));

y = [gamma_s_dot1, gamma_dot1_s, gamma_dot1_dot2, gamma_dot2_dot1, gamma_dot2_d, gamma_d_dot2];
