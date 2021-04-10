function [charge_density, current] = master_equation0(Vd, Vg, n)
%master_equation0 computes the master equation 
%   [currentDensity, current] = master_equation0(Vd, Vg, n) returns current 
%   and charge_density for n electrons, drain voltage Vd and gate voltage Vg
%   charge_density is a vector of the charge density for each n
%   current is the total current calculate as the sum of charge_density


for i = -n:1:n
    %ma(N+i+1,N+i+2)  = - (gamma_L->dot + gamma_R->dot + gamma_dot->L + gamma_dot->R)
    ma(n+i+1,n+i+2) = -(f_tunnel(1,i,Vd,Vg)+f_tunnel(2,i,Vd,Vg)+f_tunnel(3,i,Vd,Vg)+f_tunnel(4,i,Vd,Vg));
    %ma(N+i+1,N+i+2) = gamma_L->dot + gamma_R->dot
    ma(n+i+1,n+i+2-1) = f_tunnel(1,i-1,Vd,Vg)+f_tunnel(4,i-1,Vd,Vg);
    %ma(N+i+1,N+i+2+1) = gamma_dot->L è gamma_dot->R
    ma(n+i+1,n+i+2+1) = f_tunnel(2,i+1,Vd,Vg)+f_tunnel(3,i+1,Vd,Vg);
end
mma = [1,zeros(1,2*n+2);...
      ma;zeros(1,2*n+2),1];
[p,o] = eig(mma);   % produces a diagonal matrix o of eigenvalues and a full
                    % matrix p whose columns are the corresponding eigenvectors  
                    % so that mma*p = p*o.
[s,pos] = min(abs(diag(o)));  % s is the minimum eigenvalue contained in 
                                % the diagonal matrix o, pos is the index
                                % of the eigenvalue 
p_n(:,1) = abs(p(:,pos))/sum(abs(p(:,pos)));    % probability that n charges are into the dot
for i = 1:length(p_n) 
    % compute the current density as the multiplication between the probability p_n
    % and the difference between tunneling rate dot->R and R->dot
    charge_density(i) = p_n(i)*(f_tunnel(1,i-n-2,Vd,Vg) - f_tunnel(2,i-n-2,Vd,Vg));
end
% the sum of currentDensity is the current
current = sum(charge_density);
end