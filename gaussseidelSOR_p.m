function [p1,k,residue] = gaussseidelSOR_p(dx,dy,p,pRHS,c1,c2,c3,c4,c5)
%% Gauss-Seidel SOR Iteration Method
p1 = p;
p_old = p1;
k = 1;

residue = 1;
tol = 10^-3;
w = 1;

while (residue>tol)
for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1
       p1(j,i)= (1-w)*p_old(j,i)+w*(-c1(j,i)*p_old(j,i+1) - c2(j,i)*p1(j,i-1) - c4(j,i)*p1(j-1,i) - c5(j,i)*p_old(j+1,i) + pRHS(j,i))/c3(j,i);
    end
end
        p1(:,1) = p1(:,2);
        p1(:,size(dx,2)) = p1(:,size(dx,2)-1);
        p1(1,:) = p1(2,:);
        p1(size(dy,2),:) = p1(size(dy,2)-1,:);

        res=zeros([size(dy,2),size(dx,2)]);

for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1
        res(j,i)= c1(j,i)*p1(j,i+1) + c2(j,i)*p1(j,i-1) + c3(j,i)*p1(j,i) + c4(j,i)*p1(j-1,i) + c5(j,i)*p1(j+1,i) - pRHS(j,i);
    end
 end

residue = mean(mean(abs(res)));


p_old = p1;
k = k+1;

end