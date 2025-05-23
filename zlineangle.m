function [angle] = zlineangle(Q)

for ii=1:length(Q)
    Qi=Q{ii};
    
    D=squareform(pdist([Qi(:,1),Qi(:,2)],'euclidean'));D(D==0)=1e9;
    
    angles=zeros(1,size(Q,1));
    for kk=1:length(Qi)
       [val ind]=min(D(kk,:));        
       angles(kk)=angol([1,0],[Qi(ind,1:2)-Qi(kk,1:2)]);
    end
    
vec=sort([angles,angles+360]);

win=ceil(length(vec)/8);


d1=diff(smooth(vec,win));
d2=smooth(diff(d1),win);
P=d2(1:end-1).*d2(2:end);
pi=find(P<=0);
d3=smooth(diff(d2),win);
% figure;plot(vec);hold on;plot(d2*1000);stem(pi,vec(pi),'r+');plot(d3*5000)

pi(vec(pi)<45)=[];
pi(vec(pi)>720-45)=[];


giradreta=vec(pi(d3(pi)>0));


angle(ii)=mean(giradreta-[1:length(giradreta)]*180);
% angle(ii)=mean(mod(giradreta,180));



end


end