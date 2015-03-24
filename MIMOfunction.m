function [ratio_js,acc_js,time_js,ratio_bp,acc_bp,time_bp,flag]=MIMOfunction(SNR,K,Nr)
%% Senario setup
a=zeros(K,1);
flag=0;


a(1,1)=rand+sqrt(-1)*rand;
a(2,1)=Nr*(rand+sqrt(-1)*rand);

Sigma=sqrt(1/(2*10^(SNR/10)));

step=1;
res=step;
res=res/180*pi;
theta=-40:step:40;% 1Hz to 200Hz with stepsize
N=length(theta);
%ind=randsample(N,K);%index of K targets
ind=[58,60];
 
off=zeros(N,1);
off_mag=step;
off(ind)=off_mag*(rand(K,1)-0.5);
%off(ind)=delta*ones(length(ind),1);
 
theta0=theta(ind)'+off(ind);%real frequency locations
lambda=20*sqrt(Sigma/2*2*log(length(theta)*2));

%% convert all the angles
theta=theta/180*pi;
off=off/180*pi;
off_mag=off_mag/180*pi;
theta0=theta0/180*pi;
step=step/180*pi;
 
%% set up the transmitters and receivers parameter
L=50;%number of samples
M=10;
c=3*10^8;% light velocity
fc=10^9;%1GHz
 
x=zeros(N,1);
x(ind)=a;
X=x*[1,off'];
 
Nt=30; %number of transmitters
Nr=10; %number of receivers
 
dt=zeros(Nt,2);%locations of transmitters
dr=zeros(Nr,2);%locations of receivers
 
for n=1:Nt
    dt(n,1)=10*(rand-0.5);
    dt(n,2)=10*(rand-0.5);
end
 
for n=1:Nr
    dr(n,1)=10*(rand-0.5);
    dr(n,2)=10*(rand-0.5);
end
 
dis_T=zeros(Nt,N);
dis_R=zeros(Nr,N);
 
for n=1:N
    for i=1:Nt
        dis_T(i,n)=dt(i,1)*cos(theta(n))+dt(i,2)*sin(theta(n));
    end
    for j=1:Nr
        dis_R(j,n)=dr(j,1)*cos(theta(n))+dr(j,2)*sin(theta(n));
    end
end
 
dist_T=zeros(Nt,K);
dist_R=zeros(Nr,K);
 
for k=1:K
    for i=1:Nt
        dist_T(i,k)=dt(i,1)*cos(theta0(k))+dt(i,2)*sin(theta0(k));
    end
    for j=1:Nr
        dist_R(j,k)=dr(j,1)*cos(theta0(k))+dr(j,2)*sin(theta0(k));
    end
end
 
W=waveform(L,Nt);%waveform need to be designed
 
%% get the received signal of all the receivers
Y=zeros(L,K);%signal received by target k at time point l
for k=1:K
    for l=1:L
        Y(l,k)=a(k)*W(l,:)*exp(sqrt(-1)*2*pi*fc/c*dist_T(:,k));
    end
end
 
Z=zeros(L,Nr);%signal received by receiver nr at time point l
 
for n=1:Nr
    for k=1:K
        Z(:,n)=Z(:,n)+exp(sqrt(-1)*2*pi*fc/c*dist_R(n,k))*Y(:,k);
    end
    Z(:,n)=Z(:,n)+Sigma*(randn(L,1)+sqrt(-1)*randn(L,1));%additive noise
end
 
phi=cell(Nr,1);
R=zeros(M,Nr);
for n=1:Nr
    D=randn(L,M)*1; %special case
    [D_tmp,noused]=qr(D);
    D=D_tmp(:,1:M);%get the tight frame
    phi{n,1}=D';%make sure that phi'*phi=I
    R(:,n)=phi{n,1}*Z(:,n);
end
 
r=R(:);
%r=r+Sigma*(randn(length(r),1)+sqrt(-1)*randn(length(r),1));
%% get the matrix formulation
V=zeros(Nt,N);
 
for n=1:N
    V(:,n)=exp(sqrt(-1)*2*pi*fc/c*dis_T(:,n));
end
 
delta=cell(Nr,1);
for i=1:Nr
    delta{i,1}=[];
    for n=1:N
        delta{i,1}=[delta{i,1},exp(sqrt(-1)*2*pi*fc/c*dis_R(i,n))*W*V(:,n)];
    end
end
 
Theta=[];
for i=1:Nr
    Theta=[Theta;phi{i,1}*delta{i,1}];
end
 
r_matrix=Theta*x;
%% get the derivative part done
diff_V=zeros(Nt,N);
 
for n=1:N
    for i=1:Nt
        diff_V(i,n)=sqrt(-1)*2*pi*fc/c*(-dt(i,1)*sin(theta(n))+dt(i,2)*cos(theta(n)))*exp(sqrt(-1)*2*pi*fc/c*dis_T(i,n));
    end
end
 
Diff=cell(Nr,N);
for i=1:Nr
    for n=1:N
        Diff{i,n}=zeros(L,N);
        Diff{i,n}(:,n)=sqrt(-1)*2*pi*fc/c*(-dr(i,1)*sin(theta(n))+dr(i,2)*cos(theta(n)))*exp(sqrt(-1)*2*pi*fc/c*dis_R(i,n))*W*V(:,n)+exp(sqrt(-1)*2*pi*fc/c*dis_R(i,n))*W*diff_V(:,n);
    end
end
 
H=cell(N+1,1);
H{1,1}=Theta;
for n=1:N
    H{n+1,1}=zeros(Nr*M,N);
    Temp=[];
    for i=1:Nr
        Temp=[Temp;phi{i,1}*Diff{i,n}];
    end
    H{n+1,1}=Temp;
end
 
A=@(x) A_fhp_rect(x,H);
At=@(y) At_fhp_rect(y,H);
 
b0=A(X(:));%reconstructed signal using 
 
%% transformation all these into real domain
r0=[real(r);imag(r)];
H{1,1}=[real(H{1,1}),-imag(H{1,1});imag(H{1,1}) real(H{1,1})];
 
for i=2:N+1
    H{i,1}=[real(H{i,1}),-imag(H{i,1});imag(H{i,1}) real(H{i,1})];
end
 
A=@(x) A_fhp_rect(x,H);
At=@(y) At_fhp_rect(y,H);
 
X_comp=[real(X);imag(X)];
b_comp=A(X_comp(:));
 
%% Using FISTA without bounded constraint
t0=cputime;
opts.maxiter=4000;
opts.acc=10^-10;
opts.print=0;
[x_fista,theta_fista]=AFISTA(r0,H,lambda,opts,step/sqrt(12));
 
x_fista=x_fista(1:N)+sqrt(-1)*(x_fista(N+1:2*N));
theta_fista=(theta_fista(1:N)+theta_fista(N+1:2*N))/2;
[x_fista,theta_fista]=merge(step/2,x_fista,theta_fista,theta,res,0.2);
 
[dummy,ind_est]=sort(abs(x_fista),'descend');
ind_wrong=ind_est(K+1:N);
ind_est=ind_est(1:K);
%ratio0=10*log10(norm(x_al(ind_est))^2/norm(x_al(ind_wrong))^2*(N-K)/K);
ratio_js=min(abs(x_fista(ind_est)))/max(abs(x_fista(ind_wrong)));
 
loc_est=theta(ind_est)'+theta_fista(ind_est);
acc0=0;
for i=1:K
    minval=pi;
    for j=1:K
        if abs(loc_est(i)-theta0(j))<minval
            minval=abs(loc_est(i)-theta0(j));
        end
    end
    acc0=acc0+minval;
end
acc0=acc0/K;
acc_js=acc0*180/pi;
time_js=cputime-t0;
if acc_js>1
    flag=1;
end

%% compare joint sparse algorithm with BPDN
t0=cputime;
[x_bp,theta_bp]=PP_BPDN(r0,H,lambda,step/2);
 
x_bp=x_bp(1:N)+sqrt(-1)*(x_bp(N+1:2*N));
theta_bp=(theta_bp(1:N)+theta_bp(N+1:2*N))/2;
[x_bp,theta_bp]=merge(step/2,x_bp,theta_bp,theta,res,0.2);
 
[dummy,ind_est]=sort(abs(x_bp),'descend');
ind_wrong=ind_est(K+1:N);
ind_est=ind_est(1:K);
%ratio0=10*log10(norm(x_al(ind_est))^2/norm(x_al(ind_wrong))^2*(N-K)/K);
ratio_bp=min(abs(x_bp(ind_est)))/max(abs(x_bp(ind_wrong)));
 
loc_est=theta(ind_est)'+theta_bp(ind_est);
acc0=0;
for i=1:K
    minval=pi;
    for j=1:K
        if abs(loc_est(i)-theta0(j))<minval
            minval=abs(loc_est(i)-theta0(j));
        end
    end
    acc0=acc0+minval;
end
acc0=acc0/K;
acc_bp=acc0*180/pi;
time_bp=cputime-t0;
if acc_bp>1
    flag=1;
end