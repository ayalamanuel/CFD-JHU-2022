%% Information of the code
% This is a Finite Difference code for the incompressible Navier-Stokes
% equations in 2D. The code is written to be used for a lid-driven cavity
% flow. Fractional step is used for time dependecy, the viscous term is
% solved using Crank-Nicolson and the Non-linear term is solved by
% Adams-Bashforth.
%
% The mesh generation is coded to support non-uniform mesh
clear;clc
for N = 64
%% Selecting the size of domain and grid size
Lx = 1;
Ly = 1;

Nx = N;
Ny = N;

dt = 0.001;
T = 60;
Nt = T/dt;
Re = 100;

%% Implementing boundary and initial conditions
[u,v,p,dx,dy] = BC_IC(Lx,Ly,Nx,Ny);
u0 = u;
v0 = v;

[a1,a2,a3,a4,a5,le,lw,ln,ls,b1,b2,b3,b4,b5] = coeff(dt,Re,dx,dy);

uRHS = zeros(size(dy,2),size(dx,2));
vRHS = zeros(size(dy,2),size(dx,2));
%% Main loop

s=1;
tol=10^-12;
check = 1;
Imaxmax = 1000000;
   
tic
for n = 0:dt:60
s

 % Step 1 of the fractional step --- Estimate the velocity u^*
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
        end 
    end
 
    [u_ast,k1,residue] = gaussseidelSOR_u_ast(dx,dy,u,uRHS,a1,a2,a3,a4,a5);
    [v_ast,k2,residue2] = gaussseidelSOR_u_ast(dx,dy,v,vRHS,a1,a2,a3,a4,a5);

    k1
    k2


  % Step 2 --- Estimate the pressure
    [c1,c2,c3,c4,c5] = press_coeff(Re,dt,a1,a2,a3,a4,a5);
    for j=2:size(dy,2)-1
         for i=2:size(dx,2)-1
             pRHS(j,i)= (((le(j,i)*u_ast(j,i)+(1-le(j,i))*u_ast(j,i+1)-lw(j,i)*u_ast(j,i)-(1-lw(j,i))*u_ast(j,i-1))/dx(i))+...
             (((ln(j,i)*v_ast(j,i)+(1-ln(j,i))*v_ast(j+1,i))-ls(j,i)*v_ast(j,i)-(1-ls(j,i))*v_ast(j-1,i))/dy(j)))/dt;
         end
    end
    [p1,k3,residue3] = gaussseidelSOR_p(dx,dy,p,pRHS,c1,c2,c3,c4,c5);

    k3
  % Step 3 --- Correct the velocity with the pressure
    u_new=zeros(size(dy,2),size(dx,2));
    u_new(:,1) = 0;
    u_new(:,size(dx,2)) = 0;
    u_new(1,:) = 1;
    u_new(size(dy,2),:) = 0;
    v_new(:,1) = 0;
    v_new(:,size(dx,2)) = 0;
    v_new(1,:) = 0;
    v_new(size(dy,2),:) = 0;
    
    for j=2:size(dy,2)-1
        for i=2:size(dx,2)-1
          u_new(j,i) = u_ast(j,i) - (dt/dx(i))*(le(j,i)*p1(j,i) + (1-le(j,i))*p1(j,i+1) - lw(j,i)*p1(j,i) - (1-lw(j,i))*p1(j,i-1));
          v_new(j,i) = v_ast(j,i) - (dt/dy(j))*(ln(j,i)*p1(j,i) + (1-ln(j,i))*p1(j+1,i) - ls(j,i)*p1(j,i) - (1-ls(j,i))*p1(j-1,i));
        end
    end

%     [Div_max,Div_mean] = Div_vel(dx,dy,u_new,v_new)

    check(s) = mean(mean(abs(u_new-u)));
    mean(mean(abs(p-p1)))
    s = s+1;

    u0 = u;
    u = u_new;
    v0 = v;
    v = v_new;
    p = p1;

end
toc
end

