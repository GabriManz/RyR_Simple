function [IDi,I0o] = clusterRyrs(Qi,DX,Im,plt)
% Clusteritza coordenades de RyRs segons zlines.
% Els que no pertanyen a cap grup sels assigna grup zero.
% Imatge
% La imatge no es necessita, serveix només per veure el resultat.
if(nargin<4),plt=0;else plt=1;end


dmin=2.5/DX;% distància per incloure veins
amax=15;%  màxima desviació en graus respecte langle mitjà dels zlines


%     Troba direccio mitja dels zlines
% I0=zeros(max(Qi(:,2)),max(Qi(:,1)));
 I0=zeros(size(Im)); 
    D=squareform(pdist([Qi(:,1),Qi(:,2)],'euclidean'));D(D==0)=1e9;
     R=length(Qi); 
    angles=zeros(1,R);
  
    for kk=1:R
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
Zang=mean(giradreta-[1:length(giradreta)]*180);
Zvec=[cos(Zang*3.1416/180) sin(Zang*3.1416/180)];

% mira quins cumpleixen les condicions inicials
links=zeros(size(D));
for kk=1:R-1
    for kkk=kk+1:R
        if(D(kk,kkk)<dmin)
            angle=angol(Zvec,[Qi(kk,1:2)-Qi(kkk,1:2)]);
            if(angle>180),angle=360-angle;end
            if(angle>90),angle=180-angle;end
            if(angle<amax)
               links(kk,kkk)=1; 
            end
        end
    end
end
 [a,b]=find(links==1);
% figure;imagesc(Inorm(:,:,ii));axis image; hold on;
% for kk=1:length(a)
%      h=line([Qi(a(kk),1) Qi(b(kk),1)],[Qi(a(kk),2) Qi(b(kk),2)]);set(h,'color','w');   
% end

% agrupa linies
% for kk=1:length(a)
%    I0(min([Qi(a(kk),2) Qi(b(kk),2)]):max([Qi(a(kk),2) Qi(b(kk),2)]),min([Qi(a(kk),1) Qi(b(kk),1)]):max([Qi(a(kk),1) Qi(b(kk),1)]))=1;
% end

for ii=1:length(a)
   
    or=Qi(a(ii),:);fi=Qi(b(ii),:);
    [allx,ally] = puntsEnmig([or(1) fi(1)],[or(2) fi(2)]);
    for jj=1:length(allx)
       
        I0(ally(jj),allx(jj))=1;
        
    end
    
end


I0=bwlabel(I0);
I0o=I0;
%I0o=zeros(size(Im));
%I0o(1:max(Qi(:,2)),1:max(Qi(:,1)))=I0;
IDi=zeros(R,1);
for kk=1:R
    IDi(kk)=I0(Qi(kk,2),Qi(kk,1));
end

if(plt==1)
figure;imagesc(Im);axis image; hold on;
cm=superjet(max(IDi)+1,'lines');


for kk=1:R
    h=plot(Qi(kk,1),Qi(kk,2),'.');set(h,'markersize',30,'color',cm(IDi(kk)+1,:));
end
end

end
