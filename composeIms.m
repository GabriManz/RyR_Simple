function [Isum,Ch] =composeIms(folder,tagRyR,fr,tn,S)
% pesca tn imatges composades cadascuna de fr imatges

namesR=[];
% Si ha escollit més d'una imatge
if isempty(S),
namesR=[dir([folder,'\*' tagRyR '*.tif']);dir([folder,'\*' tagRyR '*.png']);dir([folder,'\*' tagRyR '*.jpg'])]; % Imatges RyRs
else
    for ii=1:size(S,1)
        namesR(ii).name=S(1,ii);
    end
end

I=imread([folder,'\',namesR(1).name]);
gg=class(I);transf=0;
switch gg(5)
    case '8'% uint8
   %nofemrus     
    case '1'% unit16
    I=double(I);I=I/(2^16);I=uint8(255*I);  transf=1;  
    case 'l'% double
    warning('Unexpected file format, please mail Carme.')  
end
Ch=2;
if (max(size(size(I)))>=3),
    [v,Ch]=max([sum(sum(I(:,:,1)));sum(sum(I(:,:,2)));sum(sum(I(:,:,3)))]);
else
    % ARREGLAR!
    ind=find(tagRyR=='1');
    if ~isempty(ind), Ch=1;end
end
if fr>length(namesR), fr=length(namesR); end

Iwin=zeros(size(I,1),size(I,2),fr);
Isum=zeros(size(I,1),size(I,2),tn);

cdq=ceil((length(namesR)-fr)/tn);
rang0=[1:fr]; ranga=[];
if (tn==1)&&(fr==1), % Agafa imatge central
    I=imread([folder '\' namesR(round(length(namesR)/2)).name]);
    if transf==1,
        I=double(I);I=I/(2^16);I=uint8(255*I);
    end
    if (max(size(size(I)))>=3), I=I(:,:,Ch); end
    Isum=I;
else
    for jj=1:tn% per cada imatge que volem de output
        rang=rang0+round((jj-1)*cdq); % Agafa la primera del rang --> canviar pq agafi la del mig?
        for kak=1:length(rang)% composada de fr imatges
            kk=rang(kak);
            if(isempty(intersect(kk,ranga)))
                I=imread([folder '\' namesR(kk).name]);
                if transf==1,
                    I=double(I);I=I/(2^16);I=uint8(255*I);
                end
                if (max(size(size(I)))>=3),I=I(:,:,Ch);end
                Iwin(:,:,kak)=I;
            else
                Iwin(:,:,kak)=Iwin(:,:,ranga==kk);
            end
        end
        Isum(:,:,jj)=sum(Iwin,3)/fr;
        %figure;imagesc(Isum);axis image
        ranga=rang;
    end
end


end