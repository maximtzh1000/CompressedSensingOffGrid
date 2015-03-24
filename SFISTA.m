function [s,p,niter]=SFISTA(b,A,W,lambda,r,opts,s_init,p_init)
%  Author: Zhao Tan
% This code is desinged for the bounded joint sparse recovery when signal
% magnitude is known to be nonzero
% Inputs:
% -b  the measurement
% -A  the orginal matrix
% -W  the perturbation matrix 
% -lambda tradeoff between l2 and l1 norm
% -opts the options for this algorithm
% -opts.C number of continuation
% -opts.r changing rate of continuation
% -opts.acc the accuracy we want to obtain
% -opts.rhof the final value of rho
% -opts.maxiter the maximum iteration number
% Outputs:
% -s the reconstructed signal
% -niter number of iteratios used to achieve the desired accuracy

maxiter=opts.maxiter;
acc=opts.acc;
mu=opts.mu;
print=opts.print;
fmean = realmin/10;

n=size(A,2);
M=[A,W];
s_old=s_init;
p_old=p_init;
 
u=s_old;
v=p_old;
niter=0;
  
t_old=1;

MtM=M'*M;
Mtb=M'*b;
L=norm(M)^2+1/mu;
    
iter_num=0;
OK=0;
time=[];

    while(iter_num<maxiter)
        t=(1+sqrt(1+4*t_old^2))/2;
        
        normX=norms([u,v],2,2);
                
        uTemp=u./normX.*max(normX-lambda*mu,0);
        vTemp=v./normX.*max(normX-lambda*mu,0);
        ind=find(normX<=10^-4);
        uTemp(ind)=zeros(length(ind),1);
        vTemp(ind)=zeros(length(ind),1);
        du=1/mu*(u-uTemp);
        dv=1/mu*(v-vTemp);
        
        xTemp=[u;v]-(MtM*[u;v]-Mtb+[du;dv])/L;
        
        sTemp=xTemp(1:n);
        pTemp=xTemp(n+1:2*n);
        
        % Projection part
        ind=find(sTemp+r*pTemp<0 & sTemp-r*pTemp<0);
        pTemp(ind)=zeros(length(ind),1);
        sTemp(ind)=zeros(length(ind),1);
        
        ind=find(pTemp>r*sTemp & sTemp+r*pTemp>0);
        vecTemp=[1;r]*[1,r]*[sTemp(ind)';pTemp(ind)']/(1+r^2);
        sTemp(ind)=vecTemp(1,:)';
        pTemp(ind)=vecTemp(2,:)';
        
        ind=find(sTemp-r*pTemp>0 & pTemp<-r*sTemp);
        vecTemp=[1;-r]*[1,-r]*[sTemp(ind)';pTemp(ind)']/(1+r^2);
        sTemp(ind)=vecTemp(1,:)';
        pTemp(ind)=vecTemp(2,:)';
        
        s=sTemp;
        p=pTemp;
        
        value_L2=0.5*(M*[s;p]-b)'*(M*[s;p]-b);
        f_val=value_L2+lambda*sum(norms([s,p],2,2));

        u=s+(t_old-1)/t*(s-s_old);
        v=p+(t_old-1)/t*(p-p_old);
        
        qp = abs(f_val - mean(fmean))/mean(fmean);%stop test
        if qp <= acc && OK; break;end
        if qp <= acc && ~OK; OK=1; end
        fmean = [f_val,fmean];
 
        if (length(fmean) > 10) fmean = fmean(1:10);end      
        iter_num=iter_num+1;
        if(print)
            fprintf('iter= %5d fval=%10.10f value = %10.10f\n',iter_num,f_val,norm(s_old-s)^2+norm(p_old-p)^2);
        end
        t_old=t;
        s_old=s;
        p_old=p;
       
        time=[time,cputime];
    end
    niter=niter+iter_num;