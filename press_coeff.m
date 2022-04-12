function [c1,c2,c3,c4,c5] = press_coeff(Re,dt,a1,a2,a3,a4,a5)

c1 = -2*(Re/dt)*a1;
c2 = -2*(Re/dt)*a2;
c4 = -2*(Re/dt)*a4;
c5 = -2*(Re/dt)*a5;
c3 = -2*(Re/dt)*(a3-1);


end