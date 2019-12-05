% Matlab program source for numerical simulation of Master equation
% in single electron transistor

%1. first step, the following physical constant and device parameters are defined as follows.
clc;
clear all;
close all;
% Definition of Physical constant
q=1.602e-19;    % electronic charge (C)
kb=1.381e-23;   % Boltzman constant (J/K)
% Definition of Device parameters
c1=1.0e-20;     % tunnel capacitor C1 (F)
c2=2.1e-19;     % tunnel capacitor C2 (F)
cg=1.0e-18;     % gate capacitor Cg (F)
ctotal=c1+c2+cg;    % total capacitance (F)
mega=1000000;   % definition of mega=106
r1=15*mega;     % tunnel resistance R1 (Ohm)
r2=250*mega;    % tunnel resistance R2 (Ohm)

%2. The values of external parameters (V, V_G, Q_0 and T) is given. Here,
%the V_G , Q_0 and T are kept a constant while the V is varied from Vmin to Vmax, as follows:
Vg=0;   % gate voltage (V)
q0=0;   % background charge q0 is assumed to be zero
temp=10;    % temperature T (K)
vmin=-0.5;  % drain voltage minimum Vmin (V)
vmax=0.5;   % drain voltage maximum Vmax (V)
NV=1000;    % number of grid from Vmin to Vmax
dV=(vmax-vmin)/NV;  % drain voltage increment of each grid point
for iv=1:NV     % loop start for drain voltage
    V(iv)=vmin+iv*dV;   % drain voltage in each grid point
    % Note that loop end for drain voltage is located in the end of this program source
    
    %3. Calculation of ?F, as follows:
    Nmin=-20;    % minimum number of N (charge number in dot)
    Nmax=20;     % maximum number of N (charge number in dot)
    for ne=1:Nmax-Nmin  % loop start for N
        n=Nmin+ne;   % N charge number in dot
        %Calculation of deltaE
        dE1p=q/ctotal*(0.5*q+(n*q-q0)-(c2+cg)*V(iv)+cg*Vg);
        dE1n=q/ctotal*(0.5*q-(n*q-q0)+(c2+cg)*V(iv)-cg*Vg);
        dE2p=q/ctotal*(0.5*q-(n*q-q0)-c1*V(iv)-cg*Vg);
        dE2n=q/ctotal*(0.5*q+(n*q-q0)+c1*V(iv)+cg*Vg);
        % Noted that loop end for N is located after calculation of \Gamma
        
        %4. the values of ?F are identified and then used for the calculation of \Gamma. If ?F is negative,
        % will be calculated by equations (26a) and (26b). However, if the ?F is positive, \Gamma is set to
        % be closed to the zero (very small). Note that the value of \Gamm ais always positive. These
        % identifications are done for four conditiond of ?F.
        
        if dE1p<0
            T1p(ne)=1/(r1*q*q)*(-dE1p)/(1-exp(dE1p/(kb*temp)));
            % ? positive in equation (26a)
        else
            T1p(ne)=1e-1; % ? positive is assumed to be very small
        end
        if dE1n<0
            T1n(ne)=1/(r1*q*q)*(-dE1n)/(1-exp(dE1n/(kb*temp)));
            % ? negative in equation (26a)
        else
            T1n(ne)=1e-1; % if negative is assumed to be very small
        end
        if dE2p<0
            T2p(ne)=1/(r2*q*q)*(-dE2p)/(1-exp(dE2p/(kb*temp)));
            %  positive in equation (26b)
        else
            T2p(ne)=1e-1; % if positive is assumed to be very small
        end
        
        if dE2n<0
            T2n(ne)=1/(r2*q*q)*(-dE2n)/(1-exp(dE2n/(kb*temp)));
            % ? negative in equation (26b)
        else
            T2n(ne)=1e-1; % if negative is assumed to be very small
            
        end
    end % loop end for N
    
    % 5. The p(n) of equation (28) is calculated. For this, normalization of equation (31a) must
    % be satisfied. Here, the values of p(Nmin) and p(Nmax) is assumed to be 0.01.
    
    p(1)=0.001; % ?(Nmin) is assumed to be 0.01
    p(Nmax-Nmin)=0.001; % ?(Nmax) is assumed to be 0.01
    % 6. Normalization of p is done. Here, ? p(n) is calculated.
    sum=0; % sum=0 is initial value to calculate ?
    for ne=2:Nmax-Nmin
        p(ne)=p(ne-1)*(T2n(ne-1)+T1p(ne-1))/(T2p(ne)+T1n(ne));
        % calculation of ?(N) in equation (28)
        % The conditions below are used to avoid divergence of Matlab  calculation
       
        if p(ne)>1e250
            p(ne)=1e250;
        end
        if p(ne)<1e-250
            p(ne)=1e-250;
        end
        % ---------------------
        sum=sum+p(ne);
    end
    for ne=2:Nmax-Nmin
        p(ne)=p(ne)/sum; % Normalization in equation (31b)
    end
    
    %Finally, the current is computed as follows:
    
    sumI=0; % sumI=0 is initial condition for current calculation
    
    for ne=2:Nmax-Nmin
        sumI=sumI+p(ne)*(T2p(ne)-T2n(ne));
    end
    I(iv)=q*sumI; % I in equation (32b)
end % end of drain voltage loop
figure('Name','plot of I vs V_d','NumberTitle','off');
plot(V,I); % plot of I vs V
xlabel ('Drain voltage V_d');
ylabel ('Drain current I_d');
for iv=1:NV-1
    dIdV(iv)=(I(iv+1)-I(iv))/dV; % calculation of dI/dV
end
figure('Name','plot of dI/dV vs V_d','NumberTitle','off');
plot(V(1,1:NV-1),dIdV);
xlabel ('Drain voltage V_d');
ylabel ('dI/dV');
%Coulomb Blockade Plot : Vd vs Vg
%deltaVd = e/cg;
deltaVcbp = (q/2)/(cg + c1); %usato per il fascio di rette con slope positivo da plottare
deltaVcbn = q/(2*c2); %usato per il fascio di rette con slope negativo da plottare
Vg_min = -0.5;
Vg = (Vg_min : 0.01 : -Vg_min);
for i = 1:length(Vg)
    Vcbp = ((q/2) + Vg *cg)/(cg+c2); %Vcb positive slope
    Vcbn = q/(2*c2) - (Vg * cg)/c2; %Vcb negative slope
end
Vx = (q/2)/(2*cg);  %punti in cui le rette intersecano l'asse x
numeroFasciRette = 4;
legendString=[];
figure('Name','Coulomb Blockade','NumberTitle','off');
for i = -numeroFasciRette/2 : numeroFasciRette/2
    plot (Vg, Vcbp + i*deltaVcbp);
    %area(Vg, Vd1 + i*deltaVg, max(Vd1 + (i-1)*deltaVg)); %plots Y versus X and fills the area between 0 and Y. The values in X can be numeric, datetime, duration or categorical values. 
    hold on;
    plot (Vg, Vcbn + i*deltaVcbn);
    % Identificazione dei rombi
    %%%%%%%%%%%%%% le rette si toccano per Vg = (e/2-q0)/cg &&&&&&&&&&&&&
    %pgon = polyshape([-q/c1 -Vx q/c1 Vx],[-q/ctotal 0 q/ctotal 0]);
    %plot(pgon);
    
    %tmp=find(Vd1 < (q0 + cg * i * Vx - e/2)/(ctotal-c1));
    %patch([Vg(tmp)],[Vd1(tmp)],ones(length(tmp),1),'facecolor','r')
    
    %tmp=find(y1 < 0);
    %patch([x(tmp)],[y1(tmp)],ones(length(tmp),1),'facecolor','r')
end
grid on;
xlabel ('Gate voltage V_g');
ylabel ('Coulomb Blockade voltage V_{Cb}');
hold off;
