function [IQsum, Qall, Isn] = processIms(Isum, resolution,llindar, Ch)

IQsum=zeros(size(Isum));

for ii=1:size(Isum,3)
    Im=Isum(:,:,ii); Im=double(Im);
    [s1,s2]=size(Im);
    
    % imatge original
    Imo=Im;
    
    % Im=(Im-min(min(Im)))/(max(max(Im))-min(min(Im)));
    % Corretgeix l'histograma
    % imposant que la mediana de la imatge ha de ser el fons.
    % corretgeix per lesquerra
    [aux,idk]=max(histc(reshape(Im,1,numel(Im)),0:255));idk=idk-1;
    Im=Im-idk;Im(Im<0)=0;
    % corretgeix per la dreta
    Ghh=log10(eps+histc(reshape(Im,1,numel(Im)),0:255));
    Ghh(Ghh<0)=0;
    top=find(cumsum(Ghh)/sum(Ghh)>.99,1,'first');
    Im(Im>top)=top;
    
    % Normalitza les imatges
    Im=(Im-min(min(Im)))/(max(max(Im))-min(min(Im))); % imagesc(Im);
    ImF1=Im;
    
    % ensenya( ' Estimating cell mask.');
    % mask=cellMask4(Im);
    % mask=imdilate(mask,ones(diametre,diametre));
    
    
    % convolucio inicial
    diametre=round(.2/resolution); % .4
    if(mod(diametre,2)==0),diametre=diametre+1;end
    if(diametre<5),diametre=5;end
    
    F=gaussiana2d(diametre);
    F=F/numel(F);
    ImF=conv2(Im,F,'same');%figure;imagesc(ImF);
    
    %ImF=ImF/max(max(abs(ImF)));
    %     figure(1),
    %     subplot(211),imagesc(Im); axis image;
    %     subplot(212),imagesc(ImF);axis image;
    
    L=bwlabel(ImF>llindar);%DAD=regionprops(L,'Orientation' ,'Area','MajorAxisLength','MinorAxisLength');
    
    
    % area minima dels ryrs
    %umb=ceil(diametre/2);if(umb<3)
    umb=1;
    
    % identify clusters
    % ensenya( ' Identifying clusters.');
    Q=findClusters(L,ImF,F,umb);
    
    Qall{ii}=Q;
    
    IQ=zeros(size(Im));
    for kk=1:length(Q)
        IQ(Q(kk,2),Q(kk,1))=1;
    end
    IQsum(:,:,ii)=IQ;
    Isn(:,:,ii)=ImF1;
    
end


end