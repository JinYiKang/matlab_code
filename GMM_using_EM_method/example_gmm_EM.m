clc
clear all
close all
%��һ������
mu1=[1.1,1.1];  %��ֵ
S1=[0.2,0;0,0.6];  %Э����
data1=mvnrnd(mu1,S1,200);   %������˹�ֲ�����

%%�ڶ�������
mu2=[2.25,2.25 ];
S2=[0.6,0;0,0.2];
data2=mvnrnd(mu2,S2,200);

%������������
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
title('���ݵ�ԭʼ�ֲ�');

[mean,cov,coef,p_for_GMM] = GMM_EM(data,3);
[~,label] = max(p_for_GMM,[],2);
figure(2)
scatter(data(label==1,1),data(label==1,2),'r');
hold on
scatter(data(label==2,1),data(label==2,2),'g');
hold on
scatter(data(label==3,1),data(label==3,2),'b');
title('ģ��Ԥ������ݷ���');