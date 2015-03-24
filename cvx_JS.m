function [x_js,theta]=cvx_JS(b,H,lambda,delta)

N=size(H,1)-1;
[m,n]=size(H{1,1});
W=zeros(m,n);
for i=1:N
    W=W+H{i+1,1};
end
A=H{1,1};

%l=zeros(n,1);
cvx_begin quiet
     variable x(n)
     variable p(n)
     variable B(n,2)
     minimize (sum(square_abs(b-(A*x)-(W*p)))+2*lambda*sum(norms(B,2,2)))
     subject to
             B(:,1)==x;
             B(:,2)==p./delta;
cvx_end

x_js=x;
theta=p./x;

% H_new=zeros(m,n);
% b_new=b-H{1,1}*x;
% for i=1:N
%     H_new(:,i)=H{i+1,1}*x;
% end
% theta=H_new\b_new;