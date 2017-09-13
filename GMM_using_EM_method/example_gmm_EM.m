clc
clear all
close all
%第一类数据
mu1=[1.1,1.1];  %均值
S1=[0.2,0;0,0.6];  %协方差
data1=mvnrnd(mu1,S1,200);   %产生高斯分布数据

%%第二类数据
mu2=[2.25,2.25 ];
S2=[0.6,0;0,0.2];
data2=mvnrnd(mu2,S2,200);

%第三个类数据
mu3=[-0.25,2.25 ];
S3=[0.2,0;0,0.3];
data3=mvnrnd(mu3,S3,200);
data =[data1;data2;data3];

figure(1);
scatter(data1(:,1),data1(:,2),'r');
hold on
scatter(data2(:,1),data2(:,2),'b');
hold on
scatter(data3(:,1),data3(:,2),'g');
title('数据的原始分布');

[mean,cov,coef,p_for_GMM] = GMM_EM(data,3);
[~,label] = max(p_for_GMM,[],2);
figure(2)
scatter(data(label==1,1),data(label==1,2),'r');
hold on
scatter(data(label==2,1),data(label==2,2),'g');
hold on
scatter(data(label==3,1),data(label==3,2),'b');
title('模型预测的数据分类');