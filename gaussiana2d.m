function [B] = gaussiana2d(N)


% 
%N=9;
B=zeros(N,N);
sig1=N*.2;
sig2=N*.2;
for ii=1:N
    
for jj=1:N    
    B(ii,jj)=sqrt(exp(-((ii-ceil(N/2))^2)/sig1/sig1)*exp(-((jj-ceil(N/2))^2)/sig2/sig2));
end

    
end
% 
% B=B-min(min(B));
% B=B/max(max(B));
% B=B-mean(mean(B));
% figure;imagesc(B);


end