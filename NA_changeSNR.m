% change SNR to see how it affects the JS and BPDN for nested array
clear;close all;clc;
% randomseed = RandStream('mcg16807','Seed',2);
% RandStream.setGlobalStream(randomseed);

T=1000;
Monte=50; % how many monte carlo run for each point
SNR=-10:2:10;%Number of measurements

Acc=zeros(Monte,length(SNR));% These are for the joint sparse algorithm
Time=zeros(Monte,length(SNR));

Acc0=zeros(Monte,length(SNR));% These are for the PP_BPDN algorithm
Time0=zeros(Monte,length(SNR));

for i=1:length(SNR)
    m=1;
    while(1)
        disp(['The SNR is ',num2str(SNR(i)),', and the Monte Carlo run is ',num2str(m)]);
        [Acc(m,i),Time(m,i),Acc0(m,i),Time0(m,i)]=NAfunction(SNR(i),T);
        m=m+1;
        if (m>Monte)
            break;
        end
    end
end

%% plot the results
fntsz = 16; lwdth = 1.2; %display parameter

figure(1)
plot(SNR,mean(Acc0),'-^r','linewidth',lwdth);
hold on;
plot(SNR,mean(Acc),'-xb','linewidth',lwdth);

%ylim([0,0.05]);
xlabel('SNR (dB)','fontsize',fntsz);
ylabel('DOA Estimation Error','fontsize',fntsz);
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('P-BPDN','BJS',fntsz);
set(h_legend,'FontSize',fntsz);