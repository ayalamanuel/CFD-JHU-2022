function [a1,a2,a3,a4,a5,le,lw,ln,ls,b1,b2,b3,b4,b5] = coeff(dt,Re,dx,dy)

% Alocating the data
a1 = zeros(size(dy,2),size(dx,2));
a2 = zeros(size(dy,2),size(dx,2));
a4 = zeros(size(dy,2),size(dx,2));
a5 = zeros(size(dy,2),size(dx,2));
a3 = zeros(size(dy,2),size(dx,2));

for j = 2:size(dy,2)-1
    for i = 2:size(dx,2)-1

        a1(j,i) = -dt/(2*Re)*(2/(dx(i)*(dx(i) + dx(i+1))));
        a2(j,i) = -dt/(2*Re)*(2/(dx(i)*(dx(i) + dx(i-1))));

        a5(j,i) = -dt/(2*Re)*(2/(dy(j)*(dy(j) + dy(j+1))));
        a4(j,i) = -dt/(2*Re)*(2/(dy(j)*(dy(j) + dy(j-1))));

        a3(j,i) = 1+(dt/(2*Re))*(((2/(dx(i+1)+dx(i))+(2/(dx(i-1)+dx(i)))) ...
            *(1/dx(i)))+((2/(dy(j+1)+dy(j))+(2/(dy(j-1)+dy(j))))*(1/dy(j))));
        
    end
end

for j=2:size(dy,2)-1
    for i=2:size(dx,2)-1
        le(j,i) = dx(i+1)/(dx(i)+dx(i+1));
        lw(j,i) = dx(i-1)/(dx(i)+dx(i-1));
        ln(j,i) = dy(j+1)/(dy(j)+dy(j+1));
        ls(j,i) = dy(j-1)/(dy(j)+dy(j-1));
    end
end


b1 = -a1;
b2 = -a2;
b3 = 2-a3;
b4 = -a4;
b5 = -a5;

end

       