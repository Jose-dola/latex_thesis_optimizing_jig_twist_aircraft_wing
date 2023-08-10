function a = getGaussianPosteriorLotsOfSamplesUnidimensional(f,x_training,n,x_test,nt,y_training,s,l,sdnoise) 
% regression function (real function, function being approximated)
% n is the number of training points
% nt is the number of test points
% if prior=1 then it plots prior samples

numberofsamplefunctions = 100;
k = @(x,y,s,l) s * exp(-(0.5)*((x-y)^2)/(l*l));

%n = 5;
%m = 15;
%N = n*m;
%nt = N-n;
%xtotal = linspace(0,1,N);
%I = zeros(1,N);
%for i = 1:n      
%I(m*i-floor(m/2)) = 1;
%end
%I = logical(I);
%x = xtotal(I);
%xt = xtotal(not(I));
x = x_training;
y = y_training;
xt = x_test;

clear Kxx;
for i=1:n
  for j=1:n
%    if i==j
%      Kxx(i,j) = k(x(i),x(j),s,sd);% + sdnoise;
%    else
      Kxx(i,j) = k(x(i),x(j),s,l);
%    end
  end
end
clear Kxtxt
for i=1:nt
  for j=1:nt
      Kxtxt(i,j) = k(xt(i),xt(j),s,l);% + sdnoise;
  end
end
clear Kxtx
for i=1:nt
  for j=1:n
      Kxtx(i,j) = k(xt(i),x(j),s,l);
  end
end
clear Kxxt
for i=1:n
  for j=1:nt
      Kxxt(i,j) = k(x(i),xt(j),s,l);
  end
end
K = [ Kxx Kxxt; Kxtx Kxtxt ];

xxt = [x xt];

fsize=36; fname='times';
linewidth = 2.5;
pointssize = 250; 
pointtype = '.';

%if prior == 1
% figure(); hold on;
% c = colormap(lines);
% for i = 1:numberofsamplefunctions;
%   y = mvnrnd(zeros(n+nt,1),K);
%   xy = [xxt.' y.'];
%   xysort = sortrows(xy,1);
%   C = c( mod(i+1,length(c)) , :);
%   s = scatter(xysort(:,1),xysort(:,2),pointssize,pointtype);
%   s.MarkerEdgeColor = C;
%   s.MarkerFaceColor = C;
%   plot(xysort(:,1),xysort(:,2),'Color',C);
% end
% xlabel('x (inputs)','FontName',fname,'fontsize',fsize);
% ylabel('y(x) (outputs)','FontName',fname,'fontsize',fsize);
% ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
% hold off;
%end

%% training data
%a = 1; b = 0.5; m=0.3; c=0.03;
%b1=b; m1=-m; c1=c; b2=b; m2=m; c2=c;
%f = @(x) a*x.^2 + b1*exp((-(x-m1).^2)/c1) + b2*exp((-(x-m2).^2)/c2);
%f = @(x) 3*(sin(2*pi*x)+4*x);
%y = f(x).' + randn(n,1)*sdnoise_samples;
%% ploting samples of Gaussian posterior

%L = cholesky(Kxx);
%Kxxinv = inv(L')*inv(L);
K = Kxx + sdnoise*sdnoise*eye(n);
Kinvy = K \ y;
for i = 1:nt
  KinvKxxt(:,i) = K \ Kxxt(:,i);
end
m   = Kxtx*Kinvy;
Cov = Kxtxt - Kxtx*KinvKxxt;
%L   = real(cholesky(Cov));
figu = figure('WindowState','maximized'); hold on;
for i = 1:numberofsamplefunctions;
%  u = randn(nt,1);
%  y = m + L*u
  y_predictive = mvnrnd(m,Cov);
  s = scatter(xt',y_predictive,pointssize,pointtype);
  s.MarkerEdgeColor = 'r';
  s.MarkerFaceColor = 'r';
  plot(xt',y_predictive,'r');
end
s = scatter(x.',y,round(pointssize*1.5),'o','filled');
s.MarkerEdgeColor = 'k';
s.MarkerFaceColor = 'k';
xx = linspace(min([x xt]), max([x xt]), 1000);
plot(xt,m,'k','linewidth',5);
%plot(xx.',f(xx).','k','linewidth',5);
title_str = sprintf('100 samples from Gaussian posterior',sdnoise);
title(title_str);
xlabel('x (inputs)','FontName',fname,'fontsize',fsize);
ylabel('y(x) (outputs)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
hold off;
str = sprintf('priorPosterior/100SAMPLESPosterior.fig');%,sdnoise);
savefig(str);
str = sprintf('priorPosterior/100SAMPLESPosterior.eps');%,sdnoise);
saveas(gcf,str,'epsc');

%%L = cholesky(Kxx);
%%Kxxinv = inv(L')*inv(L);
%Kxxinv = inv(Kxx + sdnoise*sdnoise*eye(n));
%m   = Kxtx*Kxxinv*y;
%Cov = Kxtxt - Kxtx*Kxxinv*Kxxt;
%L   = real(cholesky(Cov));
%figure(); hold on;
%c = colormap(lines);
%scatter(x.',y,pointssize,'o','filled');
%plot(x.',f(x).','k','linewidth',5);
%for i = 1:numberofsamplefunctions;
%  u = randn(nt,1);
%  y = m + L*u
%%  y = mvnrnd(m,Cov);
%%  xy = [xt.' y];
%%  xysort = sortrows(xy,1);
%  C = c( mod(i+1,length(c)) , :);
%  s = scatter(xt',y,pointssize,pointtype);
%  s.MarkerEdgeColor = C;
%  s.MarkerFaceColor = C;
%  plot(xt',y,'Color',C);
%end
%title_str = sprintf('l = %0.5f',l);
%title(title_str);
%xlabel('x (inputs)','FontName',fname,'fontsize',fsize);
%ylabel('y(x) (outputs)','FontName',fname,'fontsize',fsize);
%ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
%hold off;
%

end
