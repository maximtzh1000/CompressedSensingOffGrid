function [s,beta]=alternatingZhu(b,H,lambda,maxiter,acc,xinit,ratio)

N=size(H,1)-1;
[m,n]=size(H{1,1});
W=zeros(m,n);
for i=1:N
    W=W+H{i+1,1};
end
A=H{1,1};

%x_old=zeros(n,1);
x_old=xinit;
beta=zeros(n,1);
k=0;
while(1)
    k=k+1;
    M=W*diag(x_old);
    y=b-A*x_old;
    cvx_begin quiet
     variable beta(n)
     minimize (sum(square_abs(y-M*beta))+ratio*sum(square_abs(beta)))
    cvx_end
    
    M=A+W*diag(beta);
    cvx_begin quiet
     variable x(n)
     minimize (lambda*norm(x,1)+0.5*sum(square_abs(b-(M*x))))
    cvx_end
    error=sum(abs(x_old-x));
    value=0.5*norm(b-(A+W*diag(beta))*x,2)^2+lambda*norm(x,1)+0.5*ratio*norm(beta,2)^2;
    if (error<acc)|| (k>maxiter)
        break;
    end
    x_old=x;
end

s=x;