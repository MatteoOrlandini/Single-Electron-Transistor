clear all; close all;

y = [];
n = 20;
N1 = -n:n;
N2 = -n:n;
for k = 1:6
    for i = 1:(2*n+1)
        for j = 1:(2*n+1)
            F(i,j,k) = f_tunnel_2dot(k,i-n,j-n,1,0);
        end
    end
end

figure; surf(N1,N2,F(:,:,1),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,2),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,3),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,4),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,5),'FaceAlpha',0.5);
hold on; surf(N1,N2,F(:,:,6),'FaceAlpha',0.5);

