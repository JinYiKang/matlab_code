clc
close all
clear all
x = 0:0.01:2;
y = sin(2*pi*x) + (2*rand(1,size(x,2)) - 1)*0.2;

[w,alpha,beta,K] = autoLinearRegr_var(x',y');
figure;
scatter(x,y,'r');
hold on
plot(x,sin(2*pi*x),'b');
hold on
fx = zeros(length(x),K);
fx(:,1) = 1;
for i = 1:K
    fx(:,i+1) = x'.^i;
end
plot(x,fx*w,'g');