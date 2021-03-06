% change noise level and see how it affects sampling reconstruction
clear;close all;clc;
randomseed = RandStream('mcg16807','Seed',2);
RandStream.setGlobalStream(randomseed);

SNR=-10:2:10;
K=2;% how many targets
Monte=50; % how many monte carlo run for each point
Ratio=zeros(Monte,length(SNR));
Time=zeros(Monte,length(SNR));
Acc=zeros(Monte,length(SNR));% These are for the joint sparse algorithm
Ratio0=zeros(Monte,length(SNR));
Time0=zeros(Monte,length(SNR));
Acc0=zeros(Monte,length(SNR));% These are for the P_BPDN algorithm
% Ratio1=zeros(Monte,length(SNR));
% Acc1=zeros(Monte,length(SNR));% These are for the alternating minimization

for i=1:length(SNR)
    m=1;
    while(1)
        disp(['The SNR is ',num2str(SNR(i)),', and the Monte Carlo run is ',num2str(m)]);
        [Ratio(m,i),Acc(m,i),Time(m,i),Ratio0(m,i),Acc0(m,i),Time0(m,i),flag]=MIMOfunction(SNR(i),K,1);
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
plot(SNR,mean(Ratio0),'-or','linewidth',lwdth);
hold on
plot(SNR,mean(Ratio),'-xb','linewidth',lwdth);
ylim([0,80]);
xlabel('SNR','fontsize',fntsz);
ylabel('R','fontsize',fntsz);
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('P-BPDN','JS',fntsz);
set(h_legend,'FontSize',fntsz);

figure(2)
plot(SNR,mean(Acc0),'-or','linewidth',lwdth);
hold on;
plot(SNR,mean(Acc),'-xb','linewidth',lwdth);
%ylim([0,0.05]);
xlabel('SNR (dB)','fontsize',fntsz);
ylabel('DOA Estimation Accuracy','fontsize',fntsz);
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('P-BPDN','JS',fntsz);
set(h_legend,'FontSize',fntsz);