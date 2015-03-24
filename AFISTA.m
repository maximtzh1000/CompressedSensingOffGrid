function [s,theta]=AFISTA(b,H,lambda,opts,delta)
% This is a fast implementation to solve joint sparsity problem with base
% mismatch
% Inputs:
%        - b measurement vector
%        - H all the system matrix
%        - lambda the trade-off parameter
%        - opts all the parameters
% Outputs:
%        - s original signal
%        -theta grid mismatch

maxiter=opts.maxiter;
acc=opts.acc;
print=opts.print;
fmean = realmin/10;

N=size(H,1)-1;
[m,n]=size(H{1,1});
W=zeros(m,n);
for i=1:N
    W=W+H{i+1,1};
end
A=H{1,1};

P=[A,W*delta];%system matrix
L=norm(P)^2;

t_old=1;
iter_num=0;
OK=0;
x_old=P'*b;
v=x_old;

while (iter_num<maxiter)
    t=(1+sqrt(1+4*t_old^2))/2;
    x=v-P'*(P*v-b)/L; 
    x1=x(1:n); x2=x(n+1:2*n);
    normX=norms([x1,x2],2,2);
    x1=x1./normX.*max(normX-lambda/L,0);
    x2=x2./normX.*max(normX-lambda/L,0);
    x=[x1;x2];
    v=x+(t_old-1)/t*(x-x_old);
    
    f_val=0.5*(P*x-b)'*(P*x-b)+lambda*sum(norms([x1,x2],2,2));
    qp = abs(f_val - mean(fmean))/mean(fmean);%stop test
    if qp <= acc && OK; break;end
    if qp <= acc && ~OK; OK=1; end
    fmean = [f_val,fmean];
 
    if (length(fmean) > 10) fmean = fmean(1:10);end      
    iter_num=iter_num+1;
    if(print)
       fprintf('iter= %5d value = %10.10f %10.10f\n',iter_num,f_val,norm(x_old-x));
     end
     t_old=t;
     x_old=x;
end
s=x(1:n);
% theta=delta*x(n+1:2*n)./s;

H_new=zeros(m,n);
b_new=b-H{1,1}*s;
for i=1:N
    H_new(:,i)=H{i+1,1}*s;
end
theta=H_new\b_new;
