% We generate a Nested array to test the results

clear;clc;close all;
randomseed = RandStream('mcg16807','Seed',4);
RandStream.setGlobalStream(randomseed);

%% Nested Array Parameters set up
d=1/2; %distance between different transceivers
Level=25; % control the number of targets
SNR=10;
fc=17;

%Nested Array
N1=10;
N2=12;
L1=(1:N1)'*d;
L2=(1:N2)'*(N1+1)*d;

sigma_s=1;%signal power
sigma_n=sigma_s/(10^(SNR/10));%noise power
T=500;%number of time samples

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
lambda=50*sqrt(sigma_n/2*2*log(length(Range)*2));


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
y0=y0-sigma_n*f;
n=size(F,2);

bd=step/2;

%% Use Joint sparse to solve the problem
tic
% cvx_begin quiet
%      variable s(n)
%      variable p(n)
%      variable B(n,2)
%      minimize (lambda*2*sum(norms(B,2,2))+sum(square_abs(y0-(F*s)-(W*p))))
%      subject to
%              s>=0;
%              B(:,1)==s;
%              B(:,2)==p;
%              -bd*s<=p<=bd*s;
% cvx_end

opts.maxiter=20000;
opts.accf=10^-10;
opts.muf=lambda^-1*10^-8;
opts.innerstep=4000;
opts.C=5;
opts.print=0;
[s,p,dummy]=SFISTAcon(y0,F,W,lambda,bd,opts);

theta=p./s;
[s,theta]=merge(step/2,s,theta,Range,step/2,0.02);
s=s/max(s);
toc

%% find the locations of estimated targets from joint sparse Basis pursuit
[~,ind_est]=sort(abs(s),'descend');
ind_wrong=ind_est(K+1:end);
ind_est=ind_est(1:K);
ratio=min(abs(s(ind_est)))/max(abs(s(ind_wrong)));
loc_est=Range(ind_est)'+theta(ind_est);

% acc=0;
% for i=1:K
%     minval=1;
%     for j=1:K
%         if abs(loc_est(i)-loc(j))<minval
%             minval=abs(loc_est(i)-loc(j));
%         end
%     end
%     acc=acc+minval;
% end
acc=sum(abs(sort(loc_est)-sort(loc')));
acc_js=acc/K;
locEstJS=Range'+theta;
x_js=s;

%% Use BPDN to solve the problem
tic

% cvx_begin quiet
%      variable s(n)
%      variable p(n)
% %      variable a
%      minimize (2*lambda*norm(s,1)+sum(square_abs(y0-(F*s)-(W*p))))
%      subject to
%              s>=0;
% %            a>=0;
%              -bd*s<=p<=bd*s;
% cvx_end

cvx_begin quiet
     variable s(n)
     variable p(n)
     variable B(n,2)
     minimize (lambda*2*sum(norms(B,2,2))+sum(square_abs(y0-(F*s)-(W*p))))
     subject to
             s>=0;
             B(:,1)==s;
             B(:,2)==p;
             -bd*s<=p<=bd*s;
cvx_end

theta=p./s;
[s,theta]=merge(step/2,s,theta,Range,step/2,0.02);
s=s/max(s);
toc

%% find the locations of estimated targets from BPDN
[~,ind_est]=sort(abs(s),'descend');
ind_wrong=ind_est(K+1:end);
ind_est=ind_est(1:K);
ratio=min(abs(s(ind_est)))/max(abs(s(ind_wrong)));
loc_est=Range(ind_est)'+theta(ind_est);

% acc=0;
% for i=1:K
%     minval=1;
%     for j=1:K
%         if abs(loc_est(i)-loc(j))<minval
%             minval=abs(loc_est(i)-loc(j));
%         end
%     end
%     acc=acc+minval;
% end
acc=sum(abs(sort(loc_est)-sort(loc')));
acc_bp=acc/K;
locEstBP=Range'+theta;
x_bp=s;

%% plot all the results
fntsz=16;
figure(1)
subplot(2,1,1)
h=stem(loc, x,'-.r');
set(h, 'Marker', 'none');
hold on,
h=stem(locEstJS,x_js);
set(h, 'Marker', 'none');
set(gcf,'Position',[200 200 640 360]);
hold off;
xlim([-1 1]);
%xlabel('sin(\theta)','fontsize',fntsz);
title('\muBJS with continuation (Running time=4.96s)','fontsize',fntsz);

subplot(2,1,2)
h=stem(loc, x,'-.r');
set(h, 'Marker', 'none');
hold on,
h=stem(locEstBP,x_bp);
set(h, 'Marker', 'none');
hold off;
xlim([-1 1]);

xlabel('sin(\theta)','fontsize',fntsz);
title('BJS (Running time=69.03s )','fontsize',fntsz);

%% display all the accuracies
disp(['The accuracy of JS algorithm is ',num2str(acc_js)]);
disp(['The accuracy of BPDN algorithm is ',num2str(acc_bp)]);