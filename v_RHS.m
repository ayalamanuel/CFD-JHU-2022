function [vRHS] = v_RHS(dx,dy,dt,u,v,u_old,v_old,a,c,d,f,a3)
for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1

       ue(j,i) = dx(i+1)/(dx(i)+dx(i+1))*u(j,i) + dx(i)/(dx(i)+dx(i+1))*u(j,i+1);
       uw(j,i) = dx(i-1)/(dx(i)+dx(i-1))*u(j,i) + dx(i)/(dx(i)+dx(i-1))*u(j,i-1);
       un(j,i) = dy(j+1)/(dy(j)+dy(j+1))*u(j,i) + dy(j)/(dy(j)+dy(j+1))*u(j+1,i);
       us(j,i) = dy(j-1)/(dy(j)+dy(j-1))*u(j,i) + dy(j)/(dy(j)+dy(j-1))*u(j-1,i);

       ve(j,i) = dx(i+1)/(dx(i)+dx(i+1))*v(j,i) + dx(i)/(dx(i)+dx(i+1))*v(j,i+1);
       vw(j,i) = dx(i-1)/(dx(i)+dx(i-1))*v(j,i) + dx(i)/(dx(i)+dx(i-1))*v(j,i-1);
       vn(j,i) = dy(j+1)/(dy(j)+dy(j+1))*v(j,i) + dy(j)/(dy(j)+dy(j+1))*v(j+1,i);
       vs(j,i) = dy(j-1)/(dy(j)+dy(j-1))*v(j,i) + dy(j)/(dy(j)+dy(j-1))*v(j-1,i);
       

       ue_old(j,i) = dx(i+1)/(dx(i)+dx(i+1))*u_old(j,i) + dx(i)/(dx(i)+dx(i+1))*u_old(j,i+1);
       uw_old(j,i) = dx(i-1)/(dx(i)+dx(i-1))*u_old(j,i) + dx(i)/(dx(i)+dx(i-1))*u_old(j,i-1);
       un_old(j,i) = dy(j+1)/(dy(j)+dy(j+1))*u_old(j,i) + dy(j)/(dy(j)+dy(j+1))*u_old(j+1,i);
       us_old(j,i) = dy(j-1)/(dy(j)+dy(j-1))*u_old(j,i) + dy(j)/(dy(j)+dy(j-1))*u_old(j-1,i);

       ve_old(j,i) = dx(i+1)/(dx(i)+dx(i+1))*v_old(j,i) + dx(i)/(dx(i)+dx(i+1))*v_old(j,i+1);
       vw_old(j,i) = dx(i-1)/(dx(i)+dx(i-1))*v_old(j,i) + dx(i)/(dx(i)+dx(i-1))*v_old(j,i-1);
       vn_old(j,i) = dy(j+1)/(dy(j)+dy(j+1))*v_old(j,i) + dy(j)/(dy(j)+dy(j+1))*v_old(j+1,i);
       vs_old(j,i) = dy(j-1)/(dy(j)+dy(j-1))*v_old(j,i) + dy(j)/(dy(j)+dy(j-1))*v_old(j-1,i);

        a2(j,i) = -a(j,i);
        c2(j,i) = -c(j,i);
        d2(j,i) = -d(j,i);
        f2(j,i) = -f(j,i);
        a32(j,i) = 2-a3(j,i);

       % Viscous term
       vt(j,i) = a2(j,i)*v(j,i+1) + c2(j,i)*v(j,i-1) + d2(j,i)*v(j+1,i) + f2(j,i)*v(j-1,i) + a32(j,i)*v(j,i);  

       % Non-linear term
       NL(j,i) = -dt*(3/(2*dx(i))*(ve(j,i)*ue(j,i) - vw(j,i)*uw(j,i)) + 3/(2*dy(i))*((vn(j,i)^2-vs(j,i)^2)) ...
           - 1/(2*dx(i))*(ve_old(j,i)*ue_old(j,i) - vw_old(j,i)*uw_old(j,i)) - 1/(2*dy(i))*(vn_old(j,i)^2 - vs_old(j,i)^2));

       vRHS(j,i) = vt(j,i) + NL(j,i);
       
    end
end 
end

