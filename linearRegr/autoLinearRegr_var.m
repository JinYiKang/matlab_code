function [wf,alphaf,betaf,K] = autoLinearRegr_var(X,T,KMAX)
%input:
%   X:N*1 vector
%   T:N*1 vector
%   KMAX：the maximum degree of polynomial，it has a default vaule 12
%output:
%   wf:parameters of the selected model
%   alphaf:precision parameter in a prior Gaussian distribution of w.
%   betaf: precision parameter in a Gaussian distribution of T.
%   K:The selected degree of polynomial.
%This function is only for two dimension data points！
%writen by JinYiKang 2017/9/21.
if nargin == 2
    KMAX = 12;
end
W = zeros(KMAX+1,KMAX);
A = zeros(KMAX,1);
B = zeros(KMAX,1);
lb_for_k = [];
for i=1:KMAX
    [w,alpha,beta,lbmax] = LinearRegr_var(X,T,i);
    W(1:(i+1),i) = w;
    A(i) = alpha;
    B(i) = beta;
    lb_for_k = [lb_for_k,lbmax];
end
figure;
plot(1:length(lb_for_k),lb_for_k);
title('模型不同阶数的最大lower bound');
xlabel('模型阶数');
ylabel('最大lower bound');

[~,label] = max(lb_for_k,[],2);
wf = W(:,label);
wf = wf(wf ~= 0);
alphaf = A(label);
betaf = B(label);
K = label;
fprintf('选定的模型阶数为%d\n',K);

[~,~,~,~,lb_plot] = LinearRegr_var(X,T,K);
figure;
plot(2:length(lb_plot),lb_plot(2:end));
title('选定阶数模型的lower bound 随迭代的变化');
xlabel('迭代次数');
end