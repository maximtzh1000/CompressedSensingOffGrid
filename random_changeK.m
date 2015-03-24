% change sparsity to see how it affects FISTA and AM methods
clear;close all;clc;
% randomseed = RandStream('mcg16807','Seed',2);
% RandStream.setGlobalStream(randomseed);

K=2:2:12;% how many targets
Monte=50; % how many monte carlo run for each point
M=50;%Number of measurements
Sigma=0.1;%magintude of additive noise
sigma=1;

Ratio=zeros(Monte,length(K));
Acc=zeros(Monte,length(K));% These are for the joint sparse algorithm
Time=zeros(Monte,length(K));
Ratio0=zeros(Monte,length(K));% These are for the PP_BPDN algorithm
Acc0=zeros(Monte,length(K));
Time0=zeros(Monte,length(K));

for i=1:length(K)
    m=1;
    while(1)
        disp(['The Sparsity K is ',num2str(K(i)),', and the Monte Carlo run is ',num2str(m)]);
        [Ratio(m,i),Acc(m,i),Time(m,i),Ratio0(m,i),Acc0(m,i),Time0(m,i)]=randomfunction(K(i),M,Sigma,sigma);
        m=m+1;
        if (m>Monte)
            break;
        end
    end
end

%% plot the results
fntsz = 16; lwdth = 1.2; %display parameter

figure(1)
plot(K,mean(Acc0),'-xb','linewidth',lwdth);
hold on;
plot(K,mean(Acc),'-^r','linewidth',lwdth);

%ylim([0,0.05]);
xlabel('Sparsity','fontsize',fntsz);
ylabel('Mismatch Estimation Accuarcy','fontsize',fntsz);
%set(gcf,'Position',[200 200 640 360]);
h_legend=legend('AM','JS',fntsz);
set(h_legend,'FontSize',fntsz);

figure(2)
plot(K,mean(Ratio0),'-^r','linewidth',lwdth);
hold on;
plot(K,mean(Ratio),'-xb','linewidth',lwdth);

%ylim([0,0.05]);
xlabel('Sparsity','fontsize',fntsz);
ylabel('Signal Reconstruction Error','fontsize',fntsz);
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('AM','JS',fntsz);
set(h_legend,'FontSize',fntsz);