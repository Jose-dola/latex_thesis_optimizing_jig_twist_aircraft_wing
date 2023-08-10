a = 2;
c1 = 10;
c2 = 1;
v1 = [a,1]';
v2 = [-1,a]';
V = [v1 v2];
W = 1/(1+a^2) * V*[c1 0; 0 c2]*V';

n=1e3;
maxx = 1e1;
x = linspace(-maxx,maxx,n);
X = x'*ones(1,n);
Y = X';
clear Z;
for i=1:n
  for j=1:n
    v = [X(i,j) Y(i,j)]';
    Z(i,j) = sqrt(v'*W*v);
  end
end

fsize=36; fname='times'; linewidth=3;
figure()
[C,h] = contour(X,Y,Z,'DisplayName','a=2, c1=10, c2=1','ShowText','on','LineWidth',3);
clabel(C,h,'FontSize',26);
xlabel('v_1','FontName',fname,'fontsize',fsize);
ylabel('v_2','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
legend show;
