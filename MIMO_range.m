% change noise level and see how it affects sampling reconstruction
clear;close all;clc;
randomseed = RandStream('mcg16807','Seed',2);
RandStream.setGlobalStream(randomseed);

Range=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
K=2;% how many targets
SNR=10;
Monte=50; % how many monte carlo run for each point
Ratio=zeros(Monte,length(Range));
Time=zeros(Monte,length(Range));
Acc=zeros(Monte,length(Range));% These are for the joint sparse algorithm
Ratio0=zeros(Monte,length(Range));
Time0=zeros(Monte,length(Range));
Acc0=zeros(Monte,length(Range));% These are for the P_BPDN algorithm
% Ratio1=zeros(Monte,length(SNR));
% Acc1=zeros(Monte,length(SNR));% These are for the alternating minimization

for i=1:length(Range)
    m=1;
    while(1)
        disp(['The Dnaymic Range is ',num2str(Range(i)),', and the Monte Carlo run is ',num2str(m)]);
        [Ratio(m,i),Acc(m,i),Time(m,i),Ratio0(m,i),Acc0(m,i),Time0(m,i),flag]=MIMOfunction(SNR,K,Range(i));
        m=m+1;
        if (m>Monte)
            break;
        end
        if(flag==1)
            m=m-1;
        end
    end
end

%% plot the results
fntsz = 16; lwdth = 1.2; %display parameter

figure(1)
plot(Range,mean(Ratio0),'-or','linewidth',lwdth);
hold on
plot(Range,mean(Ratio),'-xb','linewidth',lwdth);
ylim([0,80]);
xlabel('Dynamic Range','fontsize',fntsz);
ylabel('R','fontsize',fntsz);
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('P-BPDN','JS',fntsz);
set(h_legend,'FontSize',fntsz);

figure(2)
plot(Range,mean(Acc0),'-or','linewidth',lwdth);
hold on;
plot(Range,mean(Acc),'-xb','linewidth',lwdth);
%ylim([0,0.05]);
xlabel('Dynamic Range','fontsize',fntsz);
ylabel('DOA Estimation Accuracy','fontsize',fntsz);
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('P-BPDN','JS',fntsz);
set(h_legend,'FontSize',fntsz);