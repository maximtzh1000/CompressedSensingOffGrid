function [x_js,theta]=cvx_BJS(b,H,lambda,delta,bd)

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
     variable y(n)
     variable p(n)
     variable B(n,2)
     minimize (sum(square_abs(b-(A*x)-(W*p)+(A*y)))+2*lambda*sum(norms(B,2,2)))
     subject to
             l<=x;
             l<=y;
             B(:,1)==x-y;
             B(:,2)==p./delta;
             -bd*(x+y)<=p<=bd*(x+y)
             
cvx_end
x=x-y;
x_js=x;
theta=p./x;
