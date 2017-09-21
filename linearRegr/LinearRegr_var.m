function [w,alpha,beta,lbmax,lb_plot] = LinearRegr_var(X,T,KMAX,toPlot)
%input:
%   X:N*1 vector
%   T:N*1 vector
%   KMAX：the degree of polynomial
%   toPlot:whether to draw a plot of lower bound.
%output:
%   w:linear model parameter
%   alpha：precision parameter in a prior Gaussian distribution of w.
%   beta: precision parameter in a Gaussian distribution of T.
%   lbmax:max vaule of lower bound.
%   lb_plot:values of lower bound during the iteration.
%This function is only for two dimension data points!
%writen by JinYiKang 2017/9/21.
if nargin < 3
    fprintf('wrong number of input!');
else if nargin < 4
        toPlot = 0;
    end
end

[N,~] = size(X);
%process with X
fx = zeros(N,KMAX+1);
fx(:,1) = 1;
for i = 1:KMAX
    fx(:,i+1) = X.^i;
end
%prior distribtion parameters 
a0 = 1e-10;b0 = 1e-10;
c0 = 1e-10;d0 = 1e-10;
%initialize parameters
a = a0;b = b0;
c = c0;d = d0;
lb = 0;lb_plot = [];
thresh = 0.01;
MAXITER = 200;
while MAXITER > 0
    lb_plot = [lb_plot,lb];
    Ealpha = a/b;
    Ebeta = c/d;
    iSn = Ealpha*eye(KMAX+1,KMAX+1) + Ebeta*fx'*fx;
    Sn = inv(iSn);
    Mn = Ebeta*Sn*fx'*T;
    a = a0 + (KMAX+1)/2;
    btmp = Mn'*Mn + trace(Sn);
    b = b0 + btmp;
    c = c0 + N/2;
    dtmp = (T - fx*Mn)'*(T - fx*Mn) + trace(fx'*fx*Sn);
    d = d0 + 0.5*dtmp;
    
    %lower bound
    v1 = N/2*(psi(c) - log(d) - log(2*pi)) - c/d/2*dtmp;
    v2 = -(KMAX+1)/2*log(2*pi) + (KMAX+1)/2*(psi(a) - log(b)) -a/b/2*btmp;
    v3 = a0*log(b0) + (a0-1)*(psi(a) - log(b)) - b0*a/b - log(gamma(a));
    v4 = (c0-1)*(psi(c) - log(d)) - d0*c/d + c0*log(d) - log(gamma(c0));
    
    v5 = 0.5*log(det(Sn)) + (KMAX+1)/2*(1+log(2*pi));
    v6 = log(gamma(a)) - (a-1)*psi(a) - log(b) + a;
    v7 = (c-1)*psi(c) + log(d) - c - log(gamma(c));
    
    if abs(v1 + v2 + v3 + v4 - v5 -v6 - v7 - lb) > thresh
        lb = v1 + v2 + v3 + v4 - v5 -v6 - v7;
    else
        break;
    end
    MAXITER = MAXITER - 1;
end
if toPlot
    figure;
    plot(2:length(lb_plot),lb_plot(2:end));
    xlabel('迭代次数');
    title('lower bound变化');
end
w = Mn;
alpha = Ealpha;
beta = Ebeta;
lbmax = lb;
end