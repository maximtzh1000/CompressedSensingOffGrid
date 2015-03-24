function [x_js,theta_js]=merge(bd,x_js,theta_js,grid,res,r)
%this function is used to merge the ambuigity in the DOA estimation problem
%Inputs:  -bd the bound for mismatch
%         -x_js reconstructed signal
%         -theta_js reconstructed mismatch
%         -grid orginal grid
%         -res resolution
%         -r to what level we care about
%Outputs: -x_js merged reconstructed signal
%         -theta_js merged reconstructed mismatch

%% the first part is dealing with ambiguity from the reconstruction 
N=length(x_js);

for n=1:N
    if (theta_js(n)<-bd)
        theta_js(n)=-bd;
    elseif (theta_js(n)>bd)
        theta_js(n)=bd;
    end
end

for n=1:N-1
    if theta_js(n)>0 && theta_js(n+1)<0 %condition for merge
       theta=(abs(x_js(n))*(theta_js(n)-bd)+abs(x_js(n+1))*(theta_js(n+1)+bd))/(abs(x_js(n))+abs(x_js(n+1)));
       mag=sqrt(abs(x_js(n))^2+abs(x_js(n+1))^2);
       %mag=x_js(n)+x_js(n+1);
       if theta<0
           theta_js(n)=bd+theta;
           x_js(n)=mag;
           theta_js(n+1)=0;
           x_js(n+1)=0;
       else
           theta_js(n+1)=-bd+theta;
           x_js(n+1)=mag;
           theta_js(n)=0;
           x_js(n)=0;
       end
    end
end


%% the second part is dealing with the resolution which can minimize the coherence between the dictionary
while(1)
    [val]=max(abs(x_js));
    ind=find(abs(x_js)>r*val); %find the indices that we need to focus on
    loc_est=grid(ind)'+theta_js(ind);
    L=length(ind);
    flag=0;
    for i=1:L
        for j=1:i
            if i~=j
                if abs(loc_est(i)-loc_est(j))<res
                    theta=(loc_est(i)*abs(x_js(ind(i)))+loc_est(j)*abs(x_js(ind(j))))/(abs(x_js(ind(i)))+abs(x_js(ind(j))));
                    mag=abs(x_js(ind(i)))+abs(x_js(ind(j)));
                    x_js(ind(i))=0;
                    theta_js(ind(i))=0;
                    x_js(ind(j))=0;
                    theta_js(ind(j))=0;
                    loc=round((theta-grid(1))/(2*bd))+1;
                    x_js(loc)=mag;
                    theta_js(loc)=theta-grid(loc);
                    flag=1;
                end
            end
        end
    end
    if flag==0;
        break;
    end
end

