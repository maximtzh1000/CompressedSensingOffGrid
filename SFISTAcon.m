function [s,p,niter]=SFISTAcon(b,A,W,lambda,r,opts)
%  Created on Feb 1st, 2013
%  Author: Zhao Tan
%  This code implement the continuation strategy for SFISTA

optsin.print=opts.print;

n=size(A,2);
M=[A,W];
temp=M'*b;
s_ref=temp(1:n);
p_ref=temp(n+1:2*n);

mu0 =0.9*max(norms([s_ref,p_ref],2,2));
%mu0=0.01;
Gamma = (opts.muf/mu0)^(1/opts.C);
optsin.mu= mu0;
Gammat= (opts.accf/0.1)^(1/opts.C);
optsin.acc = 0.1;

maxiter=opts.maxiter;%max iteration number overall
optsin.maxiter=opts.innerstep;%max iteration number for one continuation

niter=0;

for c=1:opts.C
    if c==1
         [s_iter,p_iter,iter_num]=SFISTA(b,A,W,lambda,r,optsin,s_ref,p_ref);
    else
         [s_iter,p_iter,iter_num]=SFISTA(b,A,W,lambda,r,optsin,s_iter,p_iter);
    end
    niter=niter+iter_num;
    optsin.acc=optsin.acc*Gammat;
    optsin.mu=optsin.mu*Gamma;
    if niter>maxiter
        break;
    end
end
s=s_iter;
p=p_iter;