clear all; close all; clc;

sdnoise = 0.1;

s = 1;
sd = 0.3;
k = @(x,y,s,sd) s * exp(-(1/2)*((x-y).'*(x-y)/(sd*sd)));

n = 100;
x = linspace(0,1,n);

clear K;
for i=1:n
  for j=1:n
      K(i,j) = k(x(i),x(j),s,sd);
  end
end

numberofsamples = 4;

fsize=36; fname='times';
linewidth = 2.5;
pointssize = 600;
pointtype = '.';
figure(); hold on;
for i = 1:numberofsamples;
  y = mvnrnd(zeros(n,1),K);
  xy = [x.' y.'];
%  xysort = sortrows(xy,1);
  scatter(xy(:,1),xy(:,2),pointssize,pointtype);
end
xlabel('x (inputs)','FontName',fname,'fontsize',fsize);
ylabel('y(x) (outputs)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
hold off;

