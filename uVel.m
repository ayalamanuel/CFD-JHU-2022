figure(2)
Ghia_Re100 = csvread('uGhia_Re1000.csv');

uGhia_Re100= Ghia_Re100(end:-1:1,1);             
yGhia_Re100= Ghia_Re100(end:-1:1,2);                  

for j = 1:size(dy,2)
    y_32(j) = (j-1)*dy(j);
    y_32(size(dy,2)) = 1;
end   

hold on
if N == 16
    plot(uGhia_Re100,yGhia_Re100,'k-');
    plot(u(size(dy,2):-1:1,size(dy,2)/2),y_32,'ko')
end

if N == 32
    plot(u(size(dy,2):-1:1,size(dy,2)/2),y_32,'ks')
end

if N == 64
    plot(u(size(dy,2):-1:1,size(dy,2)/2),y_32,'k^')
end

if N == 128
    plot(uGhia_Re100,yGhia_Re100,'k-');
    plot(u(size(dy,2):-1:1,size(dy,2)/2),y_32,'ko')
end

if N == 256
    plot(u(size(dy,2):-1:1,size(dy,2)/2),y_32,'ks')
end
box on
grid on
ylabel('$y$','Interpreter','Latex','fontsize',15)
xlabel('$u$','Interpreter','Latex','fontsize',15)
ylim([0 1])
set(get(gca,'YLabel'),'rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
l = legend('$Ghia, \: et \: al$','$128 \times 128$','$256 \times 256$','Location','southwest');
set(l, 'Interpreter', 'Latex');

