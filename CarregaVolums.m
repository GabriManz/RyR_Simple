function [tagr,Ch,a,b,V,V2,nnr,RF]=CarregaVolums(folder,RF,vid,DT)
% Llegeix dades i carrega volums vídeos


tagr='c2'; Ch=1; Sr=dir([folder '\*' tagr '*.tif']);
if isempty(Sr), tagr='ch01'; Sr=dir([folder '\*' tagr '*.tif']);end
I=imread([folder,'\',Sr(1).name]);
[a,b,c]=size(I); nnr=length(Sr);
gg=class(I);transf=0;
switch gg(5)
    case '8'% uint8
        %nofemrus
    case '1'% unit16
        I=double(I);I=I/(2^16);I=uint8(255*I);  transf=1;
    case 'l'% double
        warning('Unexpected file format, please mail Carme.')
end
V=zeros([a,b,nnr]);
for ii=1:nnr
    %     [Sr(ii).name]
    I=double(imread([folder,'\',Sr(ii).name]));
    if transf==1, I=double(imread([folder,'\',Sr(ii).name])); I=I/(2^16);I=uint8(255*I); end
    if size(I,3)>1, I=I(:,:,Ch); end
    % Video:
    V(:,:,ii)=I;
end

% figure(1),for ii=1:length(Sr),imagesc(V(:,:,ii)); axis image; pause(0.1);end

% Norm video:
V2=zeros(size(V)); MM=quantile2(V(:),.99); mm=min(min(min(V)));
for ii=1:nnr, V2(:,:,ii)=(V(:,:,ii)-mm)/(MM-mm);end; V2(V2>1)=1;

% Make the movie:
if vid==1,
    ind=find(folder=='\');namefile=['video_',folder(ind(end)+1:end)];
    makemovie6(V.*.99,folder,DT,1,1,3,0, Ch, namefile, 1);
    ind=find(folder=='\'); namefile=['video_norm_',folder(ind(end)+1:end)];
    makemovie6(V2.*.99,folder,DT,1,1,5,0, Ch, namefile, 0);
end

RF=[RF '_Events']; if (exist([folder,'\',RF],'dir')==0), mkdir([folder,'\',RF]);end
save([folder,'\',RF,'\volums.mat'],'V2','V');


end