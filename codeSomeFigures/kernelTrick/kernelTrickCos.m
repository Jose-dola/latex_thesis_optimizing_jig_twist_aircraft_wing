clear all; close all; clc;

minx = -1;
maxx = 1;
n = 1000;
m = 10;

%a = 4;
%b = 1;
%c = 1;
%d = 1;
%f = @(x) a*x.^3 + b*x.^2 + c*x + d;
%x = linspace(minx,maxx,n);

a = 3;
b = 2;
c = 1;
f = @(x) a*x.*x + b*x + c;
x = linspace(minx,maxx,n);

xt = rand(1,m)*2-1;
pol  = polyfit(xt,f(xt),1); 
regLine = @(x) pol(1)*x + pol(2);
fittype = fittype('a0 + a1*cos(x-2.8)');
p = fit(xt',f(xt)',fittype,'StartPoint',[1,1])

fsize=36; fname='times';
functionlinecolor = 'b';
regressionmodellinecolor = 'm';
regresionlinecolor = 'r';
pointscolor = 'k';
pointssize = 800;
pointtype = '.';
linewidth = 2.5;
%
%figure()
%title('underfit / overfit');
%plot(x,f(x),'DisplayName','function','Color',functionlinecolor,'LineWidth',linewidth); hold on;

figure()
plot(x,f(x),'DisplayName','function: 3*x^2+2*x+1','Color',functionlinecolor,'LineWidth',linewidth);
hold on;
scatter(xt,f(xt),pointssize,pointtype,pointscolor,'DisplayName','training points');
plot(x,regLine(x),'DisplayName','linear regression: f(x) = a0 + a1*x','Color',regresionlinecolor,'LineWidth',linewidth);
plot(x,p(x),'DisplayName','cos regression: f(x) = a0 + a1*cos(x-2.8)','Color',regressionmodellinecolor,'LineWidth',linewidth);
xlabel('x','FontName',fname,'fontsize',fsize);
ylabel('f(x)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
legend show;
hold off;

%%%% proyection into feature space %%%%
xmin = min(cos(xt-2.8));
xmax = max(cos(xt-2.8));
%ymin = min(sin(x));
%ymax = max(sin(x));
rx = linspace(xmin,xmax,10);
%ry = linspace(ymin,ymax,10);                                                                            
%X = rx'*ones(1,10);
%Y = ry'*ones(1,10);
%Y = Y';
frx = p.a0 + p.a1*rx;

xGraph = @(x) cos(x-2.8);

figure()
plot(rx,frx,'Color','r','DisplayName','line: f(Phi) = a0 + a1*Phi','LineWidth',linewidth);
hold on;
%plot(xGraph(x),f(x),'DisplayName','line: f(x) = 3*x^2+2*x+1','Color',functionlinecolor,'LineWidth',linewidth);
%plot(xGraph(x),p(x),'DisplayName','line: f(x) = a0 + a1*cos(x-2.8)','Color',regressionmodellinecolor,'LineWidth',linewidth);
scatter(xGraph(xt),f(xt),pointssize,pointtype,pointscolor,'DisplayName','training points');
xlabel('Phi / cos(x-2.8)','FontName',fname,'fontsize',fsize);
ylabel('f(Phi) / f(x)','FontName',fname,'fontsize',fsize);
ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
legend show;
hold off;



%minx = -1;
%maxx = 1;
%k = 26;
%%xdata = rand(1,m)*(maxx-minx)+minx;
%%m4 = round(m/4);
%%xdata = rand(1,m4)*(maxx-minx)/4+minx;
%%xdata = [xdata rand(1,m4)*(maxx-minx)/1.5+minx];
%%xdata = [xdata rand(1,m-2*m4)*(maxx-minx)+minx];
%xdata = linspace(minx,maxx,k);
%p = 2*pi/(maxx-minx);
%n = 10000;
%x = linspace(minx,maxx,n);
%
%%a0 = 1;
%%a1 = 2;
%%a2 = 1;
%%b1 = 1;
%%b2 = 2;
%%displacement = 0.15;
%%r = 0.9;
%%f = @(x) -(a0+a1*cos(r*(x+0.15)*p)+b1*sin(r*(x+0.15)*p)+a2*cos(2*r*(x+0.15)*p)+b2*sin(2*r*(x+0.15)*p)) + 7
%
%%a=1; b1=0.5; m1=-0.3; c1=0.01; b2=0.5; m2=0.3; c2=0.01;
%a = 1; b = 0.5; m=0.3; c=0.03;
%b1=b; m1=-m; c1=c; b2=b; m2=m; c2=c;
%f = @(x) a*x.^2 + b1*exp((-(x-m1).^2)/c1) + b2*exp((-(x-m2).^2)/c2);
%
%%a = 2;
%%b = 0.1;
%%c = 4;
%%d = 0.01;
%%f = @(x) a*exp(b*x)+c*exp(d*x);
%
%ydata = f(xdata);
%randomnoise = (rand(1,k)*2-1)*(max(ydata)-min(ydata))*0.5*0.3;
%ydata = ydata + randomnoise;
%
%options2 = fitoptions('sin2');
%options2.Algorithm = "Levenberg-Marquardt";
%options8 = fitoptions('sin8');
%options8.Algorithm = "Levenberg-Marquardt";
%underfit = fit(xdata',ydata','sin2',options2);
%overfit  = fit(xdata',ydata','sin8',options8);
%
%%%%%PLOTS%%%%%%
%%%%%%%%%%%%%%%%
%fsize=36; fname='times';
%functionlinecolor = 'k';
%underfitcolor = 'b';
%overfitcolor = 'r';
%pointscolor = 'k';
%pointssize = 500;
%pointtype = 'x';
%linewidth = 2.5;
%
%figure()
%title('underfit / overfit');
%plot(x,f(x),'DisplayName','function','Color',functionlinecolor,'LineWidth',linewidth); hold on;
%plot(x,underfit(x),'DisplayName','underfitting','Color',underfitcolor,'LineWidth',linewidth); 
%plot(x,overfit(x),'DisplayName','overfitting','Color',overfitcolor,'LineWidth',linewidth);
%scatter(xdata,ydata,pointssize,pointtype,pointscolor,'DisplayName','training points');
%xlabel('x','FontName',fname,'fontsize',fsize);
%ylabel('f(x)','FontName',fname,'fontsize',fsize);
%ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
%legend show;
%hold off;
%
%figure()
%title('underfit');
%plot(x,f(x),'DisplayName','function','Color',functionlinecolor,'LineWidth',linewidth); hold on;
%plot(x,underfit(x),'DisplayName','underfitting','Color',underfitcolor,'LineWidth',linewidth); 
%scatter(xdata,ydata,pointssize,pointtype,pointscolor,'DisplayName','training points');
%xlabel('x','FontName',fname,'fontsize',fsize);
%ylabel('f(x)','FontName',fname,'fontsize',fsize);
%ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
%legend show;
%hold off;
%
%figure()
%title('underfit / overfit');
%plot(x,f(x),'DisplayName','function','Color',functionlinecolor,'LineWidth',linewidth); hold on;
%plot(x,overfit(x),'DisplayName','overfitting','Color',overfitcolor,'LineWidth',linewidth);
%scatter(xdata,ydata,pointssize,pointtype,pointscolor,'DisplayName','training points');
%xlabel('x','FontName',fname,'fontsize',fsize);
%ylabel('f(x)','FontName',fname,'fontsize',fsize);
%ha=gca;set(ha,'linewidth',linewidth,'FontName',fname,'FontSize',fsize,'Box','off');
%legend show;
%hold off;
%
