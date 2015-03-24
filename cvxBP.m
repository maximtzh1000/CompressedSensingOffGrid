function [s]=cvxBP(b,H,lambda,tau)

N=size(H,1);
N=N-1;
M=size(H{1,1},1);
L=size(H{1,1},2);

cvx_begin
    variable B(L, N+1);
    variable HXt(M,N+1);
    minimize (tau*norm_nuc(B)+lambda*sum(norms(B,2,2)))
    subject to
        b==sum(HXt,2);
        HXt(:,1)==H{1,1}*B(:,1);
        for i=2:N+1
            HXt(:,i)==H{i,1}*B(:,i);
        end
%        Temp=sqrt(sum(B.^2,2));
cvx_end

s=B(:);