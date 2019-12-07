% f_tunnel_2dot calculates tunneling rate for a two dot transistor
%   j is j-th gamma
%   n1 is the number of electrons stored in first dot
%   n2 is the number of electrons stored in second dot
%   Vd is drain voltage
%   Vg is gate voltage

function y = f_tunnel_2dot(j,n1,n2,Vd,Vg)

%definizione delle costanti 
%e = 1.602*10^(-19);    % electronic charge (C)
e = 1;
%atto = 10^(-18);    
%Cg = 1 * atto;          %gate capacitor (F)
%Cd = 1 * atto;       %drain capacitor (F)
%Cs = 1 * atto;        %source capacitor (F)
%Cdot = 1 * atto;        %dot capacitor (dot1-dot2) (F)
Cg = 1 ;          %gate capacitor (F)
Cd = 1 ;       %drain capacitor (F)
Cs = 1 ;        %source capacitor (F)
Cdot = 1;        %dot capacitor (dot1-dot2) (F)
Ctot = Cg+Cd+Cs+Cdot; %total capacitance (F)
temp = 10; % temperature
kb = 1.381*10^(-15);   % Boltzman constant (J/K)
%mega = 10^6;     
%r = 25*mega; % tunnel resistance R1 (Ohm)
r = 1000;
Vg1 = Vg;
Vg2 = Vg;
Cg1 = Cg;
Cg2 = Cg;


%formule per due gate e due capacità di gate
%%%% Vdot1 in funzione di Vdot2
%V1 = (-n1*e*Ctot+Cg*Vg*Ctot+Cdot*(-n2*e+Cg*Vg)+Cdot*Cd*Vd)/(Ctot^2-Cdot^2); %tensione dot 1
%%%% Vdot2 in funzione di Vdot1
%V2 = (-n2*e+Cdot*V1+Cd*Vd+Cg*Vg)/Ctot; %tensione dot 2
%syms n1 e Vg1 Cg1 Vdot2 Cdot Cs n2 Vdot1 Vg2 Cg2 Vd Cd
%eqn1 = Vdot1 == (-n1*e + Vg1*Cg1 + Vdot2*Cdot)/(Cs + Cg1 + Cdot); %tensione dot 1
%eqn2 = Vdot2 ==(-n2*e + Vdot1*Cdot + Vg2*Cg2 + Vd*Cd)/(Cdot + Cg2 + Cd) ; %tensione dot 2
%y = solve(eqn1,eqn2,Vdot1,Vdot2); 
Vdot1 = (Cd*Cdot*Vd + Cd*Cg1*Vg1 + Cdot*Cg1*Vg1 + Cdot*Cg2*Vg2 + Cg1*Cg2*Vg1 - Cd*e*n1 - Cdot*e*n1 - Cdot*e*n2 - Cg2*e*n1)/(Cd*Cdot + Cd*Cg1 + Cdot*Cg1 + Cdot*Cg2 + Cg1*Cg2 + Cd*Cs + Cdot*Cs + Cg2*Cs);
Vdot2 = (Cd*Cdot*Vd + Cd*Cg1*Vd + Cdot*Cg1*Vg1 + Cdot*Cg2*Vg2 + Cg1*Cg2*Vg2 + Cd*Cs*Vd + Cg2*Cs*Vg2 - Cdot*e*n1 - Cdot*e*n2 - Cg1*e*n2 - Cs*e*n2)/(Cd*Cdot + Cd*Cg1 + Cdot*Cg1 + Cdot*Cg2 + Cg1*Cg2 + Cd*Cs + Cdot*Cs + Cg2*Cs);


%formule per un gate e una capacità di gate
%%%% Vdot1 in funzione di Vdot2
%V1 = (-n1*e*Ctot+Cg*Vg*Ctot+Cdot*(-n2*e+Cg*Vg)+Cdot*Cd*Vd)/(Ctot^2-Cdot^2); %tensione dot 1
%%%% Vdot2 in funzione di Vdot1
%V2 = (-n2*e+Cdot*V1+Cd*Vd+Cg*Vg)/Ctot; %tensione dot 2
%syms n1 e Vg1 Vdot2 Cdot Cs n2 Vdot1 Vg Vd Cd
%eqn1 = Vdot1 == (-n1*e + Vg*Cg + Vdot2*Cdot)/(Cs + Cg + Cdot); %tensione dot 1
%eqn2 = Vdot2 ==(-n2*e + Vdot1*Cdot + Vg*Cg + Vd*Cd)/(Cdot + Cg + Cd) ; %tensione dot 2
%y = solve(eqn1,eqn2,Vdot1,Vdot2); 
Vdot1 = (Vg - e*n1 + Cd*Vg + 2*Cdot*Vg + Cd*Cdot*Vd - Cd*e*n1 - Cdot*e*n1 - Cdot*e*n2)/(Cd + 2*Cdot + Cs + Cd*Cdot + Cd*Cs + Cdot*Cs + 1);
Vdot2 = (Vg - e*n2 + Cd*Vd + 2*Cdot*Vg + Cs*Vg + Cd*Cdot*Vd + Cd*Cs*Vd - Cdot*e*n1 - Cdot*e*n2 - Cs*e*n2)/(Cd + 2*Cdot + Cs + Cd*Cdot + Cd*Cs + Cdot*Cs + 1);


%notazione: dE_inizio_fine
dE_s_dot1 = -e*Vdot1+(e^2)/2*Ctot;
dE_dot1_s = e*Vdot1+(e^2)/2*Ctot;
dE_dot1_dot2 = e*Vdot1-e*Vdot2+2*(e^2)/2*Ctot;
dE_dot2_dot1 = -e*Vdot1+e*Vdot2+2*(e^2)/2*Ctot;
dE_dot2_d = -e*Vd+e*Vdot2+(e^2)/2*Ctot;
dE_d_dot2 = e*Vd-e*Vdot2+(e^2)/2*Ctot;

dE = [dE_s_dot1, dE_dot1_s, dE_dot1_dot2, dE_dot2_dot1, dE_dot2_d, dE_d_dot2];

%provvisorio 
%y = [0.001*dE(1)/(exp(5*dE(1)-1)),0.001*dE(2)/(exp(5*dE(2)-1)),0.001*dE(3)/(exp(5*dE(3)-1)),0.001*dE(4)/(exp(5*dE(4)-1)),0.001*dE(5)/(exp(5*dE(5)-1)),0.001*dE(6)/(exp(5*dE(6)-1))];
%
y = 1/(r*e^2)*(-dE(j))/(1-exp(dE(j)/(kb*temp)));
