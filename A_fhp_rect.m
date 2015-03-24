function [y]=A_fhp_rect(x,H)

n=size(H{1,1},2);
k=length(H);
X=reshape(x,n,k);

m=size(H{1,1},1);
y=zeros(m,1);

for i=1:k
    y=y+H{i,1}*X(:,i);
end