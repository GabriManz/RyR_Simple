function [mov,tn]=movRyRsvideo(DX, folder,tn)

% mirar si hi ha moviment
tagt='ch02';Ch=2;
St=dir([folder '\*' tagt '.tif']);

I=double(imread([folder '\' St(1).name]));
Vt=zeros([size(I), size(St, 1)]);
for ii=1:size(St, 1)
    I=double(imread([folder '\' St(ii).name]));
    Vt(:,:,ii)=I;
end

% figure;for ii=1:size(St,1), imagesc(Vt(:,:,ii),[0 max(max(max(Vt)))]); axis image;pause(0.05);end

% cross correlation
y=zeros(1,size(St, 1)-1); x=y;
for ii=1:size(St, 1)-1
    C = normxcorr2(Vt(:,:,ii),Vt(:,:,ii+1));
    % [max_cc, imax] = max(abs(C(:)));
    % [ypeak, xpeak] = ind2sub(size(C),imax(1));
    [ypeak, xpeak] = find(C==max(C(:)));
    y(ii)=ypeak-size(I,1);
    x(ii)=xpeak-size(I,2);
end;

% figure,
% subplot(211);plot(y);
% subplot(212);plot(x);

if max(x)<0.5/DX&&max(y)<0.5/DX,
    mov=0;
else
    tn=10*2; % One image for 0.5 s
end

ensenya(['Cell movement = ' num2str(mov)]);

% %  CLUSTER RyRs detectats diff frames ??
% if mov==1,
%     IQ=zeros([size(G.IQsum(:,:,1)),tn]);
%     for ii=1:tn
%         I=zeros([size(G.IQsum(:,:,1))]);
%         for jj=1:length(G.Q{ii})
%             Qi=G.Q{ii};
%             I(Qi(jj,3),Qi(jj,2))=1;
%         end
%         IQ(:,:,ii)=I;
%     end
%     r=ceil(.2/DX);
%     Isum=sum(IQ,3);% imagesc(Isum);
%     Isumo=Isum;Isum(Isum>0)=1;
%
%     % Hierarchical clust?
%     R=ceil(.2/DX);
%     Isum=sum(IQ,3);% imagesc(Isumo);
%     Isumo=Isum;Isum(Isum>0)=1;
%     [r,c]=find(Isum==1); pos=[r,c];
%
%     Z=linkage(pos,'ward','Euclidean'); % dendrogram(Z);
%
%     cc=cluster(Z,'cutoff',R,'criterion','distance');
%     I=zeros(size(Isumo));
%     for ii=1:length(r),I(r(ii),c(ii))=cc(ii);end
%     Qv=zeros(size(IQ));
%     for ii=1:tn, Qv(:,:,ii)=IQ(:,:,ii).*I;end
%     % imagesc(I); % max(unique(cc))
%     % figure;clf; for ii=1:tn, imagesc(IQ(:,:,ii).*I); axis image; pause(.2);end
%
%     % Agegir RyRs que ha aceptat en altres frames en cada img -->
%     % recalcular els parametres?
%     uv=unique(Qv); uv(uv==0)=[];
%     for ii=1:tn
%     u=unique(Qv(:,:,ii)); u(u==0)=[];
%     memb=ismember(uv,u); ind=[];ind=find(memb==0);
%     falt=uv(ind);
%
%     % En un for per acumular les Q:
%     as=0;
%     [Qi, IZ, pealing]=mesures(Inorm, Qi, ceil(roiR/DX), DX,as,IQsum, pg);
%     Q{ii}=Qi;
%     end
% end

end
