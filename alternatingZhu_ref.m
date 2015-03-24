function [s,beta,E_s,E_theta]=alternatingZhu_ref(b,H,lambda,maxiter,acc,xinit,ratio,x_ref,theta,ind)

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
    x_old=x;
    E_s(k)=norm(x-x_ref)/norm(x_ref);
    E_theta(k)=norm(beta(ind)-theta(ind))/norm(beta(ind));
    if (error<acc)|| (k>maxiter)
        break;
    end
end

s=x;