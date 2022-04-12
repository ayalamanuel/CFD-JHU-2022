function [pRHS] = p_RHS(dx,dy,dt,u,v)

for i=2:size(dx,2)-1
    le(i)=dx(i+1)/(dx(i)+dx(i+1));
    lw(i)=dx(i-1)/(dx(i)+dx(i-1));
end

for j=2:size(dy,2)-1
    ln(j)=dy(j+1)/(dy(j)+dy(j+1));
    ls(j)=dy(j-1)/(dy(i)+dy(j-1));
end


for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1

%        ue(j,i) = (u(j,i)*dx(i+1))/(dx(i)+dx(i+1)) + (u(j,i+1)*dx(i))/(dx(i)+dx(i+1));
%        uw(j,i) = (u(j,i)*dx(i-1))/(dx(i)+dx(i-1)) + (u(j,i-1)*dx(i))/(dx(i)+dx(i-1));
% 
%        ue(j,i) = u(j,i)*le(i) + u(j,i+1)*(1-le(i));
%        uw(j,i) = u(j,i)*lw(i) + u(j,i-1)*(1-lw(i));
% 
%        vn(j,i) = dy(j+1)/(dy(j)+dy(j+1))*v(j,i) + dy(j)/(dy(j)+dy(j+1))*v(j+1,i);
%        vs(j,i) = dy(j-1)/(dy(j)+dy(j-1))*v(j,i) + dy(j)/(dy(j)+dy(j-1))*v(j-1,i);
% 
%        EW(j,i)  = (ue(j,i) - uw(j,i))/dx(i);
%        NS(j,i)  = (vn(j,i) - vs(j,i))/dy(j);
% 
%        pRHS(j,i) =  ((ue(j,i) - uw(j,i))/dx(i) + (vn(j,i) - vs(j,i))/dy(j))/dt;

       pRHS(j,i) = ((u(j,i)*le(i) + u(j,i+1)*(1-le(i)) - u(j,i)*lw(i) - u(j,i-1)*(1-lw(i)))/dx(i) + ...
                    (v(j,i)*ln(j) + v(j+1,i)*(1-ln(j)) - v(j,i)*ls(j) - v(j-1,i)*(1-ls(j)))/dy(i))/dt;

%        R3(j,i)= (((le(i)*u(j,i)+(1-le(i))*u(j,i+1)-lw(i)*u(j,i)-(1-lw(i))*u(j,i-1))/dx(i))+...
%              (((ln(j)*v(j,i)+(1-ln(j))*v(j+1,i))-ls(j)*v(j,i)-(1-ls(j))*v(j-1,i))/dy(j)))/dt;
       

             

    end
end 
end