function [x_epdn,theta]=PP_BPDN(b,H,lambda,bd)
% P_BPDN algorithm with separate positive part and negative part

N=size(H,1)-1;
[m,n]=size(H{1,1});
W=zeros(m,n);
for i=1:N
    W=W+H{i+1,1};
end
A=H{1,1};

l=zeros(n,1);
f=ones(n,1);
cvx_begin quiet
     variable x(n)
     variable p(n)
     variable y(n)
     minimize (sum(square_abs(b-(A*x)-(W*p)+(A*y)))+2*lambda*f'*(x+y))
     subject to
             l<=x;
             l<=y
             -bd*(x+y)<=p<=bd*(x+y)
             
cvx_end
x=x-y;
x_epdn=x;
theta=p./x;