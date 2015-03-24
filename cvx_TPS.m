function [x_tps,theta]=cvx_TPS(b,H,lambda,delta)

N=size(H,1)-1;
[m,n]=size(H{1,1});
W=zeros(m,n);
for i=1:N
    W=W+H{i+1,1};
end
W=W*delta;
A=H{1,1};
cvx_begin quiet
     variable x(n)
     variable p(n)
     minimize (sum(square_abs(b-(A*x)-(W*p)))+(2*lambda*norm(x,1)+2*lambda*norm(p,1)))
cvx_end
x_tps=x;
theta=delta*p./x;