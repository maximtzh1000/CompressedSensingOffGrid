function [x_epdn,theta]=cvx_BPDN(b,H,lambda,delta)

N=size(H,1)-1;
[m,n]=size(H{1,1});
W=zeros(m,n);
for i=1:N
    W=W+H{i+1,1};
end
A=H{1,1};

l=zeros(n,1);
cvx_begin quiet
     variable x(n)
     variable p(n)
     minimize (sum(square_abs(b-(A*x)-(W*p)))+(2*lambda*norm(x,1)))
     subject to
             l<=x;
             -delta*x<=p<=delta*x
             
cvx_end
x_epdn=x;
theta=p./x;