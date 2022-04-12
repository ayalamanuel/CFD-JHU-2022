figure (1)
Ghia_Re100 = csvread('vGhia_Re1000.csv');


vGhia_Re100= Ghia_Re100(:,1);             
xGhia_Re100= Ghia_Re100(:,2);                 
 

for i = 1:size(dx,2)
    x(i) = (i-1)*dx(i);
    x(size(dx,2)) = 1;
end 


hold on
if N == 16
    plot(vGhia_Re100,xGhia_Re100,'k-');
    plot(-v(size(dy,2)/2,:),x,'ko')
end

if N == 32
    plot(-v(size(dy,2)/2,:),x,'ks')
end

if N == 64
    plot(-v(size(dy,2)/2,:),x,'k^')
end

if N == 128
    plot(vGhia_Re100,xGhia_Re100,'k-');
    plot(-v(size(dy,2)/2,:),x,'ko')
end

if N == 256
    plot(-v(size(dy,2)/2,:),x,'ks')
end
box on
grid on
ylabel('$x$','Interpreter','Latex','fontsize',15)
xlabel('$v$','Interpreter','Latex','fontsize',15)
ylim([0 1])
set(get(gca,'YLabel'),'rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
l = legend('$Ghia, \: et \: al$','$128 \times 128$','$256 \times 256$','Location','southwest');
set(l, 'Interpreter', 'Latex');