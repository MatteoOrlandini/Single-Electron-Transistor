function y=master_equation(V,Vg,N);

clear ma;
for i=-N:1:N
    ma(N+i+1,N+i+2)=-(f_tunnel(1,i,V,Vg)+f_tunnel(2,i,V,Vg)+f_tunnel(3,i,V,Vg)+f_tunnel(4,i,V,Vg));
    ma(N+i+1,N+i+2-1)=f_tunnel(1,i-1,V,Vg)+f_tunnel(4,i-1,V,Vg); ma(N+i+1,N+i+2+1)=f_tunnel(2,i+1,V,Vg)+f_tunnel(3,i+1,V,Vg);
end;
mma=[1,zeros(1,2*N+2);ma;zeros(1,2*N+2),1];
[p,o]=eig(mma);
[s,posiz]=min(abs(diag(o)));
p_n(:,1)=abs(p(:,posiz))/sum(abs(p(:,posiz)));
for i=1:length(p_n) 
    contr(i)=p_n(i)*(f_tunnel(1,i-N-2,V,Vg)-f_tunnel(2,i-N-2,V,Vg));
end; 
plot(contr); 
hold on;
y=sum(contr);

