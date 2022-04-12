function [u_ast,k,residue] = gaussseidelSOR_u_ast(dx,dy,u,S,a1,a2,a3,a4,a5)
%% Gauss-Seidel SOR Iteration Method
u_ast = u;
u_old = u_ast;

residue = 1;
tol = 10^-6;
k = 1;
w = 1;

while (residue > tol)
    for j = 2:size(dy,2)-1
        for i = 2:size(dx,2)-1
            u_ast(j,i) = (1-w)*u_old(j,i) + w*(-a1(j,i)*u_old(j,i+1) - a2(j,i)*u_ast(j,i-1) - a4(j,i)*u_ast(j-1,i) - a5(j,i)*u_old(j+1,i) + S(j,i))/a3(i,j);
        end
    end
    res=zeros([size(dy,2),size(dx,2)]); 
    for j = 2:size(dy,2)-1
        for i = 2:size(dx,2)-1
            res(j,i)= a1(j,i)*u_ast(j,i+1) + a2(j,i)*u_ast(j,i-1) + a3(j,i)*u_ast(j,i) + a4(j,i)*u_ast(j-1,i) + a5(j,i)*u_ast(j+1,i) - S(j,i);
        end
    end

residue = mean(mean(abs(res)));

u_old = u_ast;
k = k+1;
end

end