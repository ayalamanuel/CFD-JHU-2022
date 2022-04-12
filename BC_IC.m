function [u,v,p,dx,dy] = BC_IC(Lx,Ly,Nx,Ny)
%% Information of the script
% This is the function in charge of creating the mesh and setting
% the boundary and initial conditions.

%% Creating the arrays for each grid sizes
% The first and last grid sizes will be zero because we are using ghost
% cells for the boundary conditions

dx = [0 (Lx/Nx)*ones(1,Nx)  0];
dy = [0 (Ly/Ny)*ones(1,Ny)  0];

%% Initial Conditions

u = zeros(size(dy,2),size(dx,2)); 
% v = zeros(size(dy,2),size(dx,2)); 

p = zeros(size(dy,2),size(dx,2)); 

%% Boundary Conditions
% 
% % Bottom BC
% u(end,:)  = 0;
% v(end,:)  = 0;
% 
% % Left BC
% u(:,1)    = 0;
% v(:,1)    = 0;
% 
% % Right BC
% u(:,end)  = 0;
% v(:,end)  = 0;
% 
% % Top BC
% u(1,:)    = 1;
% v(1,:)    = 0;

u(:,1)=0;
u(:,size(dx,2))=0;
u(1,:)=1;
u(size(dy,2),:)=0;
v(:,1)=0;
v(:,size(dx,2))=0;
v(1,:)=0;
v(size(dy,2),:)=0;
end
