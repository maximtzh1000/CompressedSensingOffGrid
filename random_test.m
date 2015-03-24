%% experiment setup
% random test fails without considering joint sparse
clear;close all;clc;
randomseed = RandStream('mcg16807','Seed',3);
RandStream.setGlobalStream(randomseed);

k=3;
m=80; %number of measurements
n=100; %dimension of the signal 

Sigma=0.1;%magintude of additive noise
sigma=1;

theta=sigma*rand(n,1);%mismatch, have ambiuguity on the zeros terms so only focus on the error of nonzero terms

%generate sparse signal
x0=randn(n,1);
ind=randsample(n,k);
x=zeros(n,1);
x(ind)=x0(ind);

X=x*[1,theta'];

%generate measurement matrix and mismatch matrices
H=cell(n+1,1);
H{1,1}=randn(m,n);
for i=2:n+1
    H{i,1}=zeros(m,n);
    H{i,1}(:,i-1)=randn(m,1);
    %H{i,1}=randn(m,n);
end

A=@(x) A_fhp_rect(x,H);
At=@(y) At_fhp_rect(y,H);

x0=X(:);
b=A(x0)+Sigma*randn(m,1);%measurement vector

lambda = 10*Sigma*sqrt(2*log(n));
r=Sigma^2/sigma^2;

%% get the result using only cvx without phase mismatches
[xcvx]=cvx_lasso(b,H{1,1},lambda);

%% using FISTA to solve the joint sparsity convex optimization problem
t0=cputime;
opts.maxiter=5000;
opts.acc=10^-10;
opts.print=0;
[x_fista,theta_fista]=AFISTA(b,H,lambda,opts,1);

ratio=norm(x_fista-x)/norm(x);
acc=sum(abs(theta_fista(ind)-theta(ind)))/k;

if acc>10
    flag=1;
end
disp(['ratio = ',num2str(ratio),' acc= ',num2str(acc)]);
disp(' ');
cputime-t0

%% Using alternating minimization method from Zhu Hao's paper
t0=cputime;
maxiter=20;
acc_zhu=10^-8;
x_init=zeros(n,1);

[x_al,theta_al,E_x,E_theta]=alternatingZhu_ref(b,H,lambda,maxiter,acc_zhu,x_init,r,x,theta,ind);

ratio0=norm(x_al-x)/norm(x);
acc0=sum(abs(theta_al(ind)-theta(ind)))/k;

if acc0>10
    flag=1;
end
disp(['ratio0 = ',num2str(ratio0),' acc0= ',num2str(acc0)]);
disp(' ');
cputime-t0

%% plot all the results
fntsz = 16; lwdth = 1.2; %display parameter

figure(1)
stem(abs(x),'xm');
hold on;
stem(abs(xcvx));
%title('reconstruction signal without considering phase mismatch');
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('Original signal','Without mismatch',fntsz);
set(h_legend,'FontSize',fntsz);

figure(2)
stem(abs(x),'xm');
hold on;
stem(abs(x_fista));
%title('reconstructed signal using joint sparsity');
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('Original signal','Joint sparse recovery',fntsz);
set(h_legend,'FontSize',fntsz);

figure(3)
stem(abs(x),'xm');
hold on;
stem(abs(x_al));
%title('reconstructed signal using alternating minimization');
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('Original signal','Alternating minimization (500 iterations)',fntsz);
set(h_legend,'FontSize',fntsz);
