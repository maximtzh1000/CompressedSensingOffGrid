function [ratio,acc,time,ratio0,acc0,time0]=randomfunction(k,m,Sigma,sigma)
%% experiment setup
n=100; %dimension of the signal 
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

%% using FISTA to solve the joint sparsity convex optimization problem
t0=cputime;
opts.maxiter=5000;
opts.acc=10^-10;
opts.print=0;
[x_fista,theta_fista]=AFISTA(b,H,lambda,opts,1);

ratio=norm(x_fista-x)/norm(x);
acc=sum(abs(theta_fista(ind)-theta(ind)))/k;

time=cputime-t0;

%% Using alternating minimization method from Zhu Hao's paper
t0=cputime;
maxiter=10;
acc_zhu=10^-8;
x_init=zeros(n,1);

[x_al,theta_al,E_x,E_theta]=alternatingZhu_ref(b,H,lambda,maxiter,acc_zhu,x_init,r,x,theta,ind);

ratio0=norm(x_al-x)/norm(x);
acc0=sum(abs(theta_al(ind)-theta(ind)))/k;

time0=cputime-t0;