clear all; close all; clc;

f = @(x) 3*(sin(2*pi*x)+4*x);
%f = @(x) (sin(12*pi*x)+6*x);

%larray = [ 0.02    0.03    0.05    0.08    0.1 0.15];%linspace(0.02,0.1,5);
sdnoisearray = linspace(0,1,5);
sdnoisearray = sdnoisearray.*0.2;

n = 5; %number of samples
m = 15;
N = n*m;
nt = N-n; %number of training points
xtotal = linspace(0,1,N);
I = zeros(1,N);
for i = 1:n      
I(m*i-floor(m/2)) = 1;
end
I = logical(I);
xtraining = xtotal(I);
xtest = xtotal(not(I));
%sdnoise_samples = 0.3;
sdnoise_samples = 0;
ytraining = f(xtraining).' + randn(n,1)*sdnoise_samples;
sdnoise = 0.05;
l = 0.15;
%for sdnoise = sdnoisearray
%  for l = larray
    getGaussianPriorAndPosteriorSamplesUniDimensional(f,xtraining,n,xtest,nt,ytraining,1,3,1,l,sdnoise);
    getGaussianPosteriorLotsOfSamplesUnidimensional(f,xtraining,n,xtest,nt,ytraining,1,l,sdnoise);
%  end
%end 

%function a = getGaussianPriorAndPosteriorSamplesUniDimensional(x_training,n,x_test,nt,y_training,prior,numberofsamplefunctions,s,l,sdnoise) 


