function [phi1,k,residue] = gaussseidelSOR_phi(dx,dy,phi,S,a1,a2,a3,a4,a5)
%% Gauss-Seidel SOR Iteration Method
phi1 = phi;
phi_old = phi1;

residue = 1;
tol = 10^-6;
k = 1;
w = 1;

while (residue > tol)
    for j = 2:size(dy,2)-1
        for i = 2:size(dx,2)-1
            phi1(j,i) = (1-w)*phi_old(j,i) + w*(-a1(j,i)*phi_old(j,i+1) - a2(j,i)*phi1(j,i-1) - a4(j,i)*phi1(j-1,i) - a5(j,i)*phi_old(j+1,i) + S(j,i))/a3(i,j);
        end
    end
    res=zeros([size(dy,2),size(dx,2)]); 
    for j = 2:size(dy,2)-1
        for i = 2:size(dx,2)-1
            res(j,i)= a1(j,i)*phi1(j,i+1) + a2(j,i)*phi1(j,i-1) + a3(j,i)*phi1(j,i) + a4(j,i)*phi1(j-1,i) + a5(j,i)*phi1(j+1,i) - S(j,i);
        end
    end

residue = mean(mean(abs(res)));

phi_old = phi1;
k = k+1;
end

end