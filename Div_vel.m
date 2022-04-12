function [Div_max,Div_mean] = Div_vel(dx,dy,u,v)

for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1

       ue(j,i) = dx(i+1)/(dx(i)+dx(i+1))*u(j,i) + dx(i)/(dx(i)+dx(i+1))*u(j,i+1);
       uw(j,i) = dx(i-1)/(dx(i)+dx(i-1))*u(j,i) + dx(i)/(dx(i)+dx(i-1))*u(j,i-1);

       vn(j,i) = dy(j+1)/(dy(j)+dy(j+1))*v(j,i) + dy(j)/(dy(j)+dy(j+1))*v(j+1,i);
       vs(j,i) = dy(j-1)/(dy(j)+dy(j-1))*v(j,i) + dy(j)/(dy(j)+dy(j-1))*v(j-1,i);

       Div_field(j,i) = ((ue(j,i) - uw(j,i))/dx(i) + (vn(j,i) - vs(j,i))/dy(j));

       Div_max = max(max(abs(Div_field)));
       Div_mean = mean(mean(abs(Div_field)));
       
    end
end 
end