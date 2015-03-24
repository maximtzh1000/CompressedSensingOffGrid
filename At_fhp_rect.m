function [x]=At_fhp_rect(y,H)

n=size(H{1,1},2);
k=length(H);
X=zeros(n,k);

for i=1:k
    X(:,i)=H{i,1}'*y;
end
x=X(:);