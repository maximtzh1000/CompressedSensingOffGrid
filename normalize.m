function [H,R]=normalize(H)

n=size(H{1,1},2);
R=zeros(2*n,1);
for i=1:n
    R(i)=norm(H{1,1}(:,i));
    H{1,1}(:,i)=H{1,1}(:,i)/R(i);
end

for i=1:n/2
    R(i+n)=norm(H{i+1,1}(:,i));
    H{i+1,1}(:,i)=H{i+1,1}(:,i)/R(i+n);
    
    R(i+n+n/2)=norm(H{i+1,1}(:,i+n/2));
    H{i+1,1}(:,i+n/2)=H{i+1,1}(:,i+n/2)/R(i+n+n/2);
end