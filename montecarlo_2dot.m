clear all; close all;

Vd = 1;
Vg = 1;
y = [];
n = 10;
N1 = -n:n;
N2 = -n:n;
for k = 1:6 %k : gamma k-esimo
    for i = 1:(2*n+1) 
        for j = 1:(2*n+1)
            F(i,j,k) = f_tunnel_2dot(k, i-(n+1), j-(n+1), Vd, Vg);
        end
    end
end
figure ('Name','2D plot','NumberTitle','off');
plot(N1, F(n+1,:,1));
hold on;
plot(N1, F(n+1,:,2));
hold on;
plot(N1, F(n+1,:,3));
hold on;
plot(N1, F(n+1,:,4));
hold on;
plot(N1, F(n+1,:,5));
hold on;
plot(N1, F(n+1,:,6));
hold on;
xlabel('Number of electron stored in first dot (N1)');
ylabel('Tunneling rate (\Gamma)');
legend ('\Gamma_{source \rightarrow dot1}', '\Gamma_{dot1 \rightarrow source}', '\Gamma_{dot1 \rightarrow dot2}', '\Gamma_{dot2 \rightarrow dot1}', '\Gamma_{dot2 \rightarrow drain}', '\Gamma_{drain2 \rightarrow dot}');

figure ('Name','3D plot','NumberTitle','off');
surf(N1,N2,F(:,:,1),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,2),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,3),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,4),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,5),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,6),'FaceAlpha',0.5);
xlabel('Number of electron stored in first dot (N1)');
ylabel('Number of electron stored in second dot (N2)');
zlabel('Tunneling rate (\Gamma)');
colorbar;