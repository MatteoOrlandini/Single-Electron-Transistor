clear all; close all;

y = [];
n = 10;
N1 = -n:n;
N2 = -n:n;

%for k = 1:6
    for i = -n : n
        for j = -n : n
            y = [y; f_tunnel_matteo(i,j,1,0)];
        end
    end
%end
%{
figure; surf(N1,N2,F(:,:,1),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,2),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,3),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,4),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,5),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,6),'FaceAlpha',0.5);
xlabel('N1');
ylabel('N2');
zlabel('\Gamma');
%}