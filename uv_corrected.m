function [u_new,v_new] = uv_corrected(dx,dy,dt,u_ast,v_ast,p1)

u_new=zeros(size(dy,2),size(dx,2));
u_new(:,1) = 0;
u_new(:,size(dx,2)) = 0;
u_new(1,:) = 1;
u_new(size(dy,2),:) = 0;
v_new(:,1) = 0;
v_new(:,size(dx,2)) = 0;
v_new(1,:) = 0;
v_new(size(dy,2),:) = 0;


for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1
% 
%        pe(j,i) = dx(i+1)/(dx(i)+dx(i+1))*p1(j,i) + dx(i)/(dx(i)+dx(i+1))*p1(j,i+1);
%        pw(j,i) = dx(i-1)/(dx(i)+dx(i-1))*p1(j,i) + dx(i)/(dx(i)+dx(i-1))*p1(j,i-1);
% 
%        pn(j,i) = dy(j+1)/(dy(j)+dy(j+1))*p1(j,i) + dy(j)/(dy(j)+dy(j+1))*p1(j+1,i);
%        ps(j,i) = dy(j-1)/(dy(j)+dy(j-1))*p1(j,i) + dy(j)/(dy(j)+dy(j-1))*p1(j-1,i);
%         
%         u_new (j,i) = u_ast(j,i) - dt/dx(i) * (pe(j,i) - pw(j,i));
%         v_new (j,i) = v_ast(j,i) - dt/dy(j) * (pn(j,i) - ps(j,i));
       

         


    end
end 
end