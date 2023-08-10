clear all; close all; clc;

% view graphs (perspective)
caz = -41.8750;
cel =  23.1951;

numberofsamplefunctions = 3;

sdnoise = 0.3;

s = 1;
l = 0.3;
k = @(x,y,s,l) s * exp(-(0.5)*((x-y)*(x-y)')/(l*l));

f = @(x,y) sin(2*pi*x) + cos(2*pi*y);

n = 5;
p = linspace(0,1,n);
X1 = p'*ones(size(p));
X2 = X1';
nt = n*10;
p = linspace(0,1,nt);
Xt1 = p'*ones(size(p));
Xt2 = Xt1';
x = [reshape(X1,n*n,1) reshape(X2,n*n,1)];
y = f(x(:,1),x(:,2)) + sdnoise*randn(n^2,1);
xt = [reshape(Xt1,nt*nt,1) reshape(Xt2,nt*nt,1)];

clear Kxx;
for i=1:n^2
  for j=1:n^2
      Kxx(i,j) = k(x(i,:),x(j,:),s,l);
  end
end
clear Kxtxt
for i=1:nt^2
  for j=1:nt^2
      Kxtxt(i,j) = k(xt(i,:),xt(j,:),s,l);% + sdnoise;
  end
end
clear Kxtx
for i=1:nt^2
  for j=1:n^2
      Kxtx(i,j) = k(xt(i,:),x(j,:),s,l);
  end
end
clear Kxxt
for i=1:n^2
  for j=1:nt^2
      Kxxt(i,j) = k(x(i,:),xt(j,:),s,l);
  end
end
K = [ Kxx Kxxt; Kxtx Kxtxt ];

fsize=46; fname='times';
linewidth = 2.5;
pointssize = 140; 
pointtype = '.';

%% plotting regression function
  figure('WindowState','maximized');
  surf(Xt1,Xt2,f(Xt1,Xt2),'FaceAlpha',0.6); hold on;
  shading flat;
  title_str = sprintf('Regression function and training points');
  title(title_str);
  xlabel('x1 (input)','FontName',fname,'fontsize',fsize);
  ylabel('x2 (input)','FontName',fname,'fontsize',fsize);
  zlabel('y(x1,x2) (output)','FontName',fname,'fontsize',fsize);
  ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
  colorbar;
  view([caz,cel]);
  grid off;
  scatter3(x(:,1),x(:,2),y,pointssize,'filled','k');
  str = sprintf('priorPosteriorNoise2dim/Function.fig');
  savefig(str);
  str = sprintf('priorPosteriorNoise2dim/Function.eps');
  saveas(gcf,str,'epsc');

%% plotting samples from Gaussian prior
for i = 1:numberofsamplefunctions;
  y_prior = mvnrnd(zeros(nt^2,1),Kxtxt);
  Y_prior = reshape(y_prior,nt,nt);
  figure('WindowState','maximized');
  surf(Xt1,Xt2,Y_prior);
  shading flat;
  title_str = sprintf('Sample from Gaussian prior');
  title(title_str);
  xlabel('x1 (input)','FontName',fname,'fontsize',fsize);
  ylabel('x2 (input)','FontName',fname,'fontsize',fsize);
  zlabel('y(x1,x2) (output)','FontName',fname,'fontsize',fsize);
  ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
  view([caz,cel]);
  grid off;
  str = sprintf('priorPosteriorNoise2dim/SamplePrior%i.fig',i);
  savefig(str);
  str = sprintf('priorPosteriorNoise2dim/SamplePrior%i.eps',i);
  saveas(gcf,str,'epsc');
end


%% plotting samples from Gaussian posterior and mean
K = Kxx + sdnoise*sdnoise*eye(n^2);
Kinvy = K \ y;
for i = 1:nt^2
  KinvKxxt(:,i) = K \ Kxxt(:,i);
end
m   = Kxtx*Kinvy;
Cov = Kxtxt - Kxtx*KinvKxxt;
figure('WindowState','maximized');
surf(Xt1,Xt2,reshape(m,nt,nt));
shading flat;
title_str = sprintf('Mean of Gaussian posterior');
title(title_str);
xlabel('x1 (input)','FontName',fname,'fontsize',fsize);
ylabel('x2 (input)','FontName',fname,'fontsize',fsize);
zlabel('y(x1,x2) (output)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
view([caz,cel]);
grid off;
str = sprintf('priorPosteriorNoise2dim/MeanPosterior.fig');
savefig(str);
str = sprintf('priorPosteriorNoise2dim/MeanPosterior.eps');
saveas(gcf,str,'epsc');
for i = 1:numberofsamplefunctions;
  y_pred = mvnrnd(m,Cov);
  Y_pred = reshape(y_pred,nt,nt);
  figure('WindowState','maximized');
  surf(Xt1,Xt2,Y_pred);
  shading flat;
  title_str = sprintf('Sample from Gaussian posterior');
  title(title_str);
  xlabel('x1 (input)','FontName',fname,'fontsize',fsize);
  ylabel('x2 (input)','FontName',fname,'fontsize',fsize);
  zlabel('y(x1,x2) (output)','FontName',fname,'fontsize',fsize);
  ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
  view([caz,cel]);
  grid off;
  str = sprintf('priorPosteriorNoise2dim/SamplePosterior%i.fig',i);
  savefig(str);
  str = sprintf('priorPosteriorNoise2dim/SamplePosterior%i.eps',i);
  saveas(gcf,str,'epsc');
end




%%figure('WindowState','maximized');
%%for i = 1:numberofsamplefunctions;
%%   y_pred = mvnrnd(m,Cov);
%%   Y_pred = reshape(y_pred,nt,nt);
%%   surf(Xt1,Xt2,Y_pred); hold on;
%%   shading flat;
%%end
%%title_str = sprintf('100 samples from Gaussian posterior');%,sdnoise);
%%xlabel('x1 (input)','FontName',fname,'fontsize',fsize);
%%ylabel('x2 (input)','FontName',fname,'fontsize',fsize);
%%zlabel('y(x1,x2) (output)','FontName',fname,'fontsize',fsize);
%%ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
%%hold off;
%
%%  C = c( mod(i+1,length(c)) , :);
%%  s = scatter(xt',y_predictive,pointssize*3,pointtype);
%%  s.MarkerEdgeColor = C;
%%  s.MarkerFaceColor = C;
%%  plot(xt',y_predictive,'Color',C,'linewidth',3);
%s = scatter(x.',y,round(pointssize*1.5),'o','filled');
%s.MarkerEdgeColor = 'k';
%s.MarkerFaceColor = 'k';
%xx = linspace(min([x xt]), max([x xt]), 1000);
%%plot(xt,m,'b','linewidth',2);
%%plot(xx.',f(xx).','k','linewidth',5);
%%title_str = sprintf('l = %0.3f, sigma^2_{noise} = %0.3f',l,sdnoise);
%title_str = sprintf('3 samples from Gaussian posterior');%,sdnoise);
%title(title_str);
%xlabel('x (input)','FontName',fname,'fontsize',fsize);
%ylabel('y(x) (output)','FontName',fname,'fontsize',fsize);
%ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
%hold off;
%%exportgraphics(gcf,str,'Resolution',300);
%
%%%L = cholesky(Kxx);
%%%Kxxinv = inv(L')*inv(L);
%%Kxxinv = inv(Kxx + sdnoise*sdnoise*eye(n));
%%m   = Kxtx*Kxxinv*y;
%%Cov = Kxtxt - Kxtx*Kxxinv*Kxxt;
%%L   = real(cholesky(Cov));
%%figure(); hold on;
%%c = colormap(lines);
%%scatter(x.',y,pointssize,'o','filled');
%%plot(x.',f(x).','k','linewidth',5);
%%for i = 1:numberofsamplefunctions;
%%  u = randn(nt,1);
%%  y = m + L*u
%%%  y = mvnrnd(m,Cov);
%%%  xy = [xt.' y];
%%%  xysort = sortrows(xy,1);
%%  C = c( mod(i+1,length(c)) , :);
%%  s = scatter(xt',y,pointssize,pointtype);
%%  s.MarkerEdgeColor = C;
%%  s.MarkerFaceColor = C;
%%  plot(xt',y,'Color',C);
%%end
%%title_str = sprintf('l = %0.5f',l);
%%title(title_str);
%%xlabel('x (input)','FontName',fname,'fontsize',fsize);
%%ylabel('y(x) (output)','FontName',fname,'fontsize',fsize);
%%ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
%%hold off;
%%
%
%end
