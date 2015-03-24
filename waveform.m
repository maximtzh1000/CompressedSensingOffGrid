function [W]=waveform(M,N)

W=zeros(M,N);
for i=1:M
    for j=1:N
        ind=randsample(4,1);
        W(i,j)=exp(sqrt(-1)*pi/4*ind);
    end
end
