function [phi1,k,residue] = gaussseidelSOR_vort2(dx,dy,phi,vortRHS,c1,c2,c3,c4,c5)
%% Gauss-Seidel SOR Iteration Method
phi1 = phi;
p_old = phi1;
k = 1;

residue = 1;
tol = 10^-3;
w = 1;

while (residue>tol)
for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1
       phi1(j,i)= (1-w)*p_old(j,i)+w*(-c1(j,i)*p_old(j,i+1) - c2(j,i)*phi1(j,i-1) - c4(j,i)*phi1(j-1,i) - c5(j,i)*p_old(j+1,i) + vortRHS(j,i))/c3(j,i);
    end
end
        phi1(:,1) = phi1(:,2);
        phi1(:,size(dx,2)) = phi1(:,size(dx,2)-1);
        phi1(1,:) = phi1(2,:);
        phi1(size(dy,2),:) = phi1(size(dy,2)-1,:);

        res=zeros([size(dy,2),size(dx,2)]);

for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1
        res(j,i)= c1(j,i)*phi1(j,i+1) + c2(j,i)*phi1(j,i-1) + c3(j,i)*phi1(j,i) + c4(j,i)*phi1(j-1,i) + c5(j,i)*phi1(j+1,i) - vortRHS(j,i);
    end
 end

residue = mean(mean(abs(res)));


p_old = phi1;
k = k+1;

end