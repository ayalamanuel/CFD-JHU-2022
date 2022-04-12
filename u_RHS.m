function [uRHS,vRHS] = u_RHS(dx,dy,dt,u,v,u0,v0,b1,b2,b3,b4,b5,le,lw,ln,ls)

for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1

  uRHS(j,i)= (b1(j,i)*u(j,i+1)+b2(j,i)*u(j,i-1)+b3(j,i)*u(j,i)+b4(j,i)*u(j-1,i)+b5(j,i)*u(j+1,i)) ...
         - dt*((1.5/dx(i))*((le(j,i)*u(j,i)+(1-le(j,i))*u(j,i+1))^2-(lw(j,i)*u(j,i)+(1-lw(j,i))*u(j,i-1))^2)...
         -(0.5/dx(i))*((le(j,i)*u0(j,i)+(1-le(j,i))*u0(j,i+1))^2-(lw(j,i)*u0(j,i)+(1-lw(j,i))*u0(j,i-1))^2)...
         + (1.5/dy(j))*((ln(j,i)*u(j,i)+(1-ln(j,i))*u(j+1,i))*(ln(j,i)*v(j,i)+(1-ln(j,i))*v(j+1,i))-(ls(j,i)*u(j,i)+(1-ls(j,i))*u(j-1,i))*(ls(j,i)*v(j,i)+(1-ls(j,i))*v(j-1,i)))...
         - (0.5/dy(j))*((ln(j,i)*u0(j,i)+(1-ln(j,i))*u0(j+1,i))*(ln(j,i)*v0(j,i)+(1-ln(j,i))*v0(j+1,i))-(ls(j,i)*u0(j,i)+(1-ls(j,i))*u0(j-1,i))*(ls(j,i)*v0(j,i)+(1-ls(j,i))*v0(j-1,i))));  

  vRHS(j,i)= (b1(j,i)*v(j,i+1)+b2(j,i)*v(j,i-1)+b3(j,i)*v(j,i)+b4(j,i)*v(j-1,i)+b5(j,i)*v(j+1,i)) ...
         - dt*((1.5/dy(j))*((ln(j,i)*v(j,i)+(1-ln(j,i))*v(j+1,i))^2-(ls(j,i)*v(j,i)+(1-ls(j,i))*v(j-1,i))^2)...
         -(0.5/dy(j))*((ln(j,i)*v0(j,i)+(1-ln(j,i))*v0(j+1,i))^2-(ls(j,i)*v0(j,i)+(1-ls(j,i))*v0(j-1,i))^2)...
         + (1.5/dx(i))*((le(j,i)*u(j,i)+(1-le(j,i))*u(j,i+1))*(le(j,i)*v(j,i)+(1-le(j,i))*v(j,i+1))-(lw(j,i)*u(j,i)+(1-lw(j,i))*u(j,i-1))*(lw(j,i)*v(j,i)+(1-lw(j,i))*v(j,i-1)))...
         - (0.5/dx(i))*((le(j,i)*u0(j,i)+(1-le(j,i))*u0(j,i+1))*(le(j,i)*v0(j,i)+(1-le(j,i))*v0(j,i+1))-(lw(j,i)*u0(j,i)+(1-lw(j,i))*u0(j,i-1))*(lw(j,i)*v0(j,i)+(1-lw(j,i))*v0(j,i-1))));      
         
%        ue(j,i) = dx(i+1)/(dx(i)+dx(i+1))*u(j,i) + dx(i)/(dx(i)+dx(i+1))*u(j,i+1);
%        uw(j,i) = dx(i-1)/(dx(i)+dx(i-1))*u(j,i) + dx(i)/(dx(i)+dx(i-1))*u(j,i-1);
%        un(j,i) = dy(j+1)/(dy(j)+dy(j+1))*u(j,i) + dy(j)/(dy(j)+dy(j+1))*u(j+1,i);
%        us(j,i) = dy(j-1)/(dy(j)+dy(j-1))*u(j,i) + dy(j)/(dy(j)+dy(j-1))*u(j-1,i);
% 
%        ve(j,i) = dx(i+1)/(dx(i)+dx(i+1))*v(j,i) + dx(i)/(dx(i)+dx(i+1))*v(j,i+1);
%        vw(j,i) = dx(i-1)/(dx(i)+dx(i-1))*v(j,i) + dx(i)/(dx(i)+dx(i-1))*v(j,i-1);
%        vn(j,i) = dy(j+1)/(dy(j)+dy(j+1))*v(j,i) + dy(j)/(dy(j)+dy(j+1))*v(j+1,i);
%        vs(j,i) = dy(j-1)/(dy(j)+dy(j-1))*v(j,i) + dy(j)/(dy(j)+dy(j-1))*v(j-1,i);
%        
% 
%        ue_old(j,i) = dx(i+1)/(dx(i)+dx(i+1))*u_old(j,i) + dx(i)/(dx(i)+dx(i+1))*u_old(j,i+1);
%        uw_old(j,i) = dx(i-1)/(dx(i)+dx(i-1))*u_old(j,i) + dx(i)/(dx(i)+dx(i-1))*u_old(j,i-1);
%        un_old(j,i) = dy(j+1)/(dy(j)+dy(j+1))*u_old(j,i) + dy(j)/(dy(j)+dy(j+1))*u_old(j+1,i);
%        us_old(j,i) = dy(j-1)/(dy(j)+dy(j-1))*u_old(j,i) + dy(j)/(dy(j)+dy(j-1))*u_old(j-1,i);
% 
%        ve_old(j,i) = dx(i+1)/(dx(i)+dx(i+1))*v_old(j,i) + dx(i)/(dx(i)+dx(i+1))*v_old(j,i+1);
%        vw_old(j,i) = dx(i-1)/(dx(i)+dx(i-1))*v_old(j,i) + dx(i)/(dx(i)+dx(i-1))*v_old(j,i-1);
%        vn_old(j,i) = dy(j+1)/(dy(j)+dy(j+1))*v_old(j,i) + dy(j)/(dy(j)+dy(j+1))*v_old(j+1,i);
%        vs_old(j,i) = dy(j-1)/(dy(j)+dy(j-1))*v_old(j,i) + dy(j)/(dy(j)+dy(j-1))*v_old(j-1,i);


%        % Viscous term
%        vt(j,i) = a2(j,i)*u(j,i+1) + c2(j,i)*u(j,i-1) + d2(j,i)*u(j+1,i) + f2(j,i)*u(j-1,i) + a32(j,i)*u(i,j);  
% 
%        % Non-linear term
%        NL(j,i) = -dt*(3/(2*dx(i))*(ue(j,i)^2-uw(j,i)^2) + 3/(2*dy(i))*((un(j,i)*vn(j,i))^2-(us(j,i)*vs(j,i))^2) ...
%            - 1/(2*dx(i))*(ue_old(j,i)^2-uw_old(j,i)^2) - 1/(2*dy(i))*((un_old(j,i)*vn_old(j,i))^2-(us_old(j,i)*vs_old(j,i))^2));
% 
%        uRHS(j,i) = vt(j,i) + NL(j,i);

    end
end 
end
