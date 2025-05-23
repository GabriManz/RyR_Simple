function [Q]=findClusters(L,ImF,F,umbral)

% [NQ]=findClusters(L,ImF,F,diametre,umbral)
% radi=(diametre-1)/2;
C=0;%Q=[];
for ii=1:max(max(L))
    %ImL(L==ii)=numel(find(L==ii));
    
    [a,b]=find(L==ii);
    if(length(a)>umbral)
        try
        A=ImF(min(a):max(a),min(b):max(b)); % Valors en una ROI quadrada de L==ii
        B=L(min(a):max(a),min(b):max(b)); % Valors de L en la ROI
        A=A.*(B==ii); % Es queda amb els valors dins de L==ii de la img ImF
        A=(A-min(min(A)))/(max(max(A))-min(min(A))); % Normalitza
        AA=imregionalmax(conv2(A,F,'same')); % Busca el màxim convolucionant la regió amb la gaussiana RyR
        AA(1,:)=0;AA(end,:)=0;AA(:,1)=0;AA(:,end)=0; % elimina els de les vores de la regió
        AA(A==0)=0;
        [aa,bb]=find(AA==1); % Busca els index del maxim
        aa=min(a)-1+aa; % els passa al num de la imatge (tenim el pix de la roi)
        bb=min(b)-1+bb;
%         for jj=length(aa):-1:1 % suprimim pics rarus en la correlacio
%             if(ImF(aa(jj),bb(jj))<=0.16),
%                 aa(jj)=[];bb(jj)=[];
%             end
%         end
        
       
        
        
        % guarda posicions
        
        for jj=1:length(aa)
            C=C+1;
            %Ctr(min(a)-1+aa(jj),min(b)-1+bb(jj))=C;
            Q(C,1)=bb(jj);Q(C,2)=aa(jj);
            %     if(length(aa)==1),
            %   Q(C,3)=DAD(ii).MajorAxisLength;Q(C,4)=DAD(ii).MinorAxisLength;Q(C,5)=-DAD(ii).Orientation;
            %     else
            
            
            %%%% ??
           % Q(C,3)=radi;Q(C,4)=radi;Q(C,5)=0;
        end
        end
    end
    
end
Normes=sqrt((Q(:,1).^2)+(Q(:,2).^2));
[~,ordre]=sort(Normes);

Q=Q(ordre,:);


% 
%  %  filtrant els paios que queden massa aprop
%   dis=squareform(pdist(Q(:,1:2)));
%   NQ=zeros(size(Q));
%   cc=0;
%   for ii=1:size(dis,1)-1
%    
%       v=dis(ii,ii+1:end);
%       a=find(v<diametre);
%       a=a+ii;
%       if(~isempty(a))
%          grup=[ii a];
%          nig=Q(grup',1:2);
%          v2=zeros(size(nig,1),1);
%          for jj=1:size(nig,1)
%              v2(jj)=ImF(nig(jj,2),nig(jj,1));
%             % if(grup(jj)==403),disp(ii);end
%          end
%          [aux ind]=max(v2);
%          for jj=1:size(nig,1)
%              
%              if(jj~=ind)
%             dis(grup(jj),:)=9999;
%             dis(:,grup(jj))=9999;
%              end
%          end
%          cc=cc+1;%if((cc==283)||(cc==284)),disp(ii);end
%          NQ(cc,1:2)=nig(ind,:);
%       else  
%           if(isempty(find(v==9999,1)))
%            cc=cc+1;%if((cc==283)||(cc==284)),disp(ii);end
%          NQ(cc,1:2)=Q(ii,1:2);
%           end
%       end
%   end
%   
%    NQ=NQ(1:cc,:);

   
%    figure;plot(Q(:,1),-Q(:,2),'r.');hold on;
%     plot(NQ(:,1),-NQ(:,2),'.')

end