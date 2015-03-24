function [acc_js,time_js,acc_bp,time_bp]=NAfunction(SNR,T)

%% Nested Array Parameters set up
d=1/2; %distance between different transceivers
Level=15; % control the number of targets
fc=17;

%Nested Array
N1=5;
N2=6;
L1=(1:N1)'*d;
L2=(1:N2)'*(N1+1)*d;

sigma_s=1;%signal power
sigma_n=sigma_s/(10^(SNR/10));%noise power

%number of targets
K = floor(Level*fc/16);
%K=15;
% nominal spacing
tnominal = (1:K);
%spike locations
tspikes = (tnominal + rand(1,K)./16).*16./(Level*fc);% OK it is the term of (-sin(phi)/2+1/2)
if tspikes(K)>1
    tspikes(K)=tspikes(K)-1;
end
tspikes=2*d*tspikes;
loc=1-tspikes/d;%locations with sin(phi)

% K=2;
% loc=[sin(20*pi/180),sin(22*pi/180)];
x=ones(K,1)*sigma_s;

step=0.01;
Range=-1:step:1;

lambda=20*sqrt(sigma_n/2*2*log(length(Range)*2));


%% Physical Model
X=(sqrt(sigma_s/2))*randn(K,T)+(sqrt(-sigma_s/2))*randn(K,T);% signal 
W=(sqrt(sigma_n/2))*randn((N1+N2),T)+(sqrt(-sigma_n/2))*randn((N1+N2),T);% noise

A0=zeros((N1+N2),K);

for i=1:N1
    for j=1:K
        A0(i,j)=exp(sqrt(-1)*L1(i)*(2*pi*loc(j)));
    end
end

for i=1:N2
    for j=1:K
        A0(i+N1,j)=exp(sqrt(-1)*L2(i)*(2*pi*loc(j)));
    end
end

Y=A0*X+W;% the received signal

R=zeros(N1+N2,N1+N2);%covariance information
for t=1:T
    R=R+Y(:,t)*Y(:,t)';
end
R=R/T;
r=R(:);

L=[L1;L2];
diffset=zeros((N1+N2)^2,1);
for i=1:(N1+N2)
    for j=1:(N1+N2)
        diffset((i-1)*(N1+N2)+j)=L(j)-L(i);
    end
end

diffset=(round(diffset*1000)/1000);%get rid of round off error
[new,I]=unique(diffset);
y=r(I);

Temp=eye(N1+N2,N1+N2);
temp=Temp(:);
temp=temp(I);

%% CS model 
F=exp(1i*2*pi*new*Range);
W=zeros(size(F));

for n=1:length(new)
    W(n,:)=1i*2*pi*new(n)*F(n,:);
end

%test the model
%use CVX to solve the construction problem
y0=[real(y);imag(y)];
F=[real(F);imag(F)];
W=[real(W);imag(W)];
f=[real(temp);imag(temp)];
n=size(F,2);

bd=step/2;
delta=1;
y0=y0-sigma_n*f;

%% Use Joint sparse to solve the problem
t0=cputime;
cvx_begin quiet
     variable s(n)
     variable p(n)
     variable B(n,2)
     minimize (lambda*2*sum(norms(B,2,2))+sum(square_abs(y0-(F*s)-(W*p))))
     subject to
             s>=0;
             B(:,1)==s;
             B(:,2)==p./delta;
             -bd*s<=p<=bd*s;
cvx_end

theta=p./s;
[s,theta]=merge(step/2,s,theta,Range,step/2,0.02);
s=s/max(s);
time_js=cputime-t0;

%% find the locations of estimated targets from joint sparse Basis pursuit
[dummy,ind_est]=sort(abs(s),'descend');
ind_wrong=ind_est(K+1:end);
ind_est=ind_est(1:K);
ratio=min(abs(s(ind_est)))/max(abs(s(ind_wrong)));
loc_est=Range(ind_est)'+theta(ind_est);

acc=0;
for i=1:K
    minval=1;
    for j=1:K
        if abs(loc_est(i)-loc(j))<minval
            minval=abs(loc_est(i)-loc(j));
        end
    end
    acc=acc+minval;
end
acc_js=acc/K;

%% Use BPDN to solve the problem
t0=cputime;

cvx_begin quiet
     variable s(n)
     variable p(n)
     minimize (2*lambda*norm(s,1)+sum(square_abs(y0-(F*s)-(W*p))))
     subject to
             s>=0;
             -bd*s<=p<=bd*s;
cvx_end

theta=p./s;
[s,theta]=merge(step/2,s,theta,Range,step/2,0.02);
s=s/max(s);
time_bp=cputime-t0;

%% find the locations of estimated targets from BPDN
[dummy,ind_est]=sort(abs(s),'descend');
ind_wrong=ind_est(K+1:end);
ind_est=ind_est(1:K);
ratio=min(abs(s(ind_est)))/max(abs(s(ind_wrong)));
loc_est=Range(ind_est)'+theta(ind_est);

acc=0;
for i=1:K
    minval=1;
    for j=1:K
        if abs(loc_est(i)-loc(j))<minval
            minval=abs(loc_est(i)-loc(j));
        end
    end
    acc=acc+minval;
end
acc_bp=acc/K;