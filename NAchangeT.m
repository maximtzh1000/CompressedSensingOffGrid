% change T to see how it affects the JS and BPDN for nested array
clear;close all;clc;
% randomseed = RandStream('mcg16807','Seed',2);
% RandStream.setGlobalStream(randomseed);

T=[200,400,600,800,1000,2000,3000,4000,5000];
Monte=50; % how many monte carlo run for each point
SNR=0;%Number of measurements

Acc=zeros(Monte,length(T));% These are for the joint sparse algorithm
Time=zeros(Monte,length(T));

Acc0=zeros(Monte,length(T));% These are for the PP_BPDN algorithm
Time0=zeros(Monte,length(T));

for i=1:length(T)
    m=1;
    while(1)
        disp(['The Number of Time Samples is ',num2str(T(i)),', and the Monte Carlo run is ',num2str(m)]);
        [Acc(m,i),Time(m,i),Acc0(m,i),Time0(m,i)]=NAfunction(SNR,T(i));
        m=m+1;
        if (m>Monte)
            break;
        end
    end
end

%% plot the results
fntsz = 16; lwdth = 1.2; %display parameter

figure(1)
plot(T,mean(Acc0),'-^r','linewidth',lwdth);
hold on;
plot(T,mean(Acc),'-xb','linewidth',lwdth);

%ylim([0,0.05]);
xlabel('Number of Time Samples','fontsize',fntsz);
ylabel('DOA Estimation Error','fontsize',fntsz);
set(gcf,'Position',[200 200 640 360]);
h_legend=legend('P-BPDN','BJS',fntsz);
set(h_legend,'FontSize',fntsz);