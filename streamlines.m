for j = 1:size(dy,2)
    y(j) = (j-1)*dy(j);
    y(size(dy,2)) = 1;
end 

for i = 1:size(dx,2)
    x(i) = (i-1)*dx(i);
    x(size(dx,2)) = 1;
end 


[X,Y] = meshgrid(x,y);

phi = zeros(size(dy,2),size(dx,2));
vortRHS = zeros(size(dy,2),size(dx,2));

for j=2:size(dy,2)-1
    for i=2:size(dx,2)-1
dvdx(j,i) = (le(j,i)*v(j,i) + (1-le(j,i))*v(j,i+1) - lw(j,i)*v(j,i) - (1-lw(j,i))*v(j,i-1))/dx(i);
dudy(j,i) = (ln(j,i)*u(j,i) + (1-ln(j,i))*u(j+1,i) - ls(j,i)*u(j,i) - (1-ls(j,i))*u(j-1,i))/dy(j);

vortRHS(j,i) = -dvdx(j,i) + dudy(j,i);
    end
end

[phi1,k,residue] = gaussseidelSOR_phi(dx,dy,phi,vortRHS,a1,a2,a3,a4,a5);
% [phi1,k,residue] = gaussseidelSOR_vort2(dx,dy,phi,vortRHS,c1,c2,c3,c4,c5);

% 
% figure (3)
% u_mag = sqrt(u.^2 + v.^2);
% contourf(X,flip(Y),(u_mag),70);
% box on
% ylabel('$y$','Interpreter','Latex','fontsize',15)
% xlabel('$x$','Interpreter','Latex','fontsize',15)
% ylim([0 1])
% set(get(gca,'YLabel'),'rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')


figure (4)
hold on
contour(X,flip(Y),qg,50)
contour(X,flip(Y),qg,[-0.0015 -0.001 -0.0007 -0.003 -4.9e-5])
% contourf(X,flip(Y),qg,50)
box on
% colormap(gray)
ylabel('$y$','Interpreter','Latex','fontsize',15)
xlabel('$x$','Interpreter','Latex','fontsize',15)
ylim([0 1])
set(get(gca,'YLabel'),'rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')