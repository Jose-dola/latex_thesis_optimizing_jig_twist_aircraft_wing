close all; clear all;

a = 1.5;
b = 3;

n = 10*1e3;
xx = linspace(0,1,1000);
x = rand(1,n);
y = 2*rand(1,n);
j = 1; 
k = 1;
for i=1:n
  if betapdf(x(i),a,b) > y(i)
    in(j,1) = x(i);
    in(j,2) = y(i);
    j = j+1;
  else
    out(k,1) = x(i);
    out(k,2) = y(i);
    k = k+1;
  end
end

fsize=36; fname='times';
linewidth = 2.5;
figure(); hold on;
s1 = scatter(in(:,1),in(:,2),'r','.');
s1.SizeData = 200;
s2 = scatter(out(:,1),out(:,2),'b','.');
s2.SizeData = 200;
p = plot(xx,betapdf(xx,a,b));
p.LineWidth = 3;
p.Color = 'k';
xlabel('x','FontName',fname,'fontsize',fsize);
ylabel('p(x)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
hold off;

figure();
histogram(in(:,1));
xlabel('x','FontName',fname,'fontsize',fsize);
ylabel('frequency','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');



