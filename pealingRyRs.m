function [suma, Qp]=pealingRyRs(IQ,DX,pg, Qi2)

radi=ceil(2.5/DX); % 2 um de diametre % 1
se=double(getnhood(strel('disk',radi)));
mask=conv2(IQ,se,'same'); mask(mask>0)=1;
L=bwlabel(mask);
if max(unique(L))>(numel(find(IQ>0))*.2),
    radi=ceil(2/DX); % 2 um de diametre % 1
    se=double(getnhood(strel('disk',radi)));
    mask=conv2(IQ,se,'same'); mask(mask>0)=1;
end
mask=imclose(mask,se); % Tanquem forats petits


% keep largest
BLUE=mask;
BLUE=bwlabel(BLUE);
kk=regionprops(BLUE,'area');
kk=cell2mat(struct2cell(kk));
[v,ik]=max(kk);
BLUE(BLUE~=ik)=0;
BLUE(BLUE~=0)=1;

%remove holes
L=bwlabel(1-BLUE);
kk=regionprops(L,'PixelList');
for ixi=1:length(kk)
    vals=kk(ixi).PixelList;
    if(isempty([find(vals(:,1)==1); find(vals(:,2)==1);find(vals(:,2)==size(L,1));find(vals(:,1)==size(L,2))]))
        BLUE(L==ixi)=1;
    end
end
mask=BLUE;

radi=ceil((2.2-.33)/DX); % contrau la mask
se=strel('disk',radi);% se=double(se.Neighborhood);
mask=imerode(mask,se);

if all(unique(mask))==1,
    radi=ceil(1.5/DX); % 2 um de diametre % 1
    se=double(getnhood(strel('disk',radi)));
    mask=conv2(IQ,se,'same'); mask(mask>0)=1;
    L=bwlabel(mask);
    if max(unique(L))>(numel(find(IQ>0))*.2),
        radi=ceil(2/DX); % 2 um de diametre % 1
        se=double(getnhood(strel('disk',radi)));
        mask=conv2(IQ,se,'same'); mask(mask>0)=1;
    end
    mask=imclose(mask,se); % Tanquem forats petits
end

% if all(unique(mask))==1,
%     radi=ceil(.8/DX); % 2 um de diametre % 1
%     se=double(getnhood(strel('disk',radi)));
%     mask=conv2(IQ,se,'same'); mask(mask>0)=1;
%     L=bwlabel(mask);
%     if max(unique(L))>(numel(find(IQ>0))*.2),
%         radi=ceil(2/DX); % 2 um de diametre % 1
%         se=double(getnhood(strel('disk',radi)));
%         mask=conv2(IQ,se,'same'); mask(mask>0)=1;
%     end
%     mask=imclose(mask,se); % Tanquem forats petits
% end

zz=0;
if all(unique(mask))==1,
    zz=1;
    % afegeix zeros
    [m,n]=size(IQ);
    IQ2=zeros(m,n+2*round(2/DX));
    IQ2(1:m,round(2/DX)+1:n+round(2/DX))=IQ;
    
    radi=ceil(.8/DX); % 2 um de diametre % 1
    se=double(getnhood(strel('disk',radi)));
    mask=conv2(IQ2,se,'same'); mask(mask>0)=1;
    L=bwlabel(mask);
    if max(unique(L))>(numel(find(IQ2>0))*.2),
        radi=ceil(1/DX); % 2 um de diametre % 1
        se=double(getnhood(strel('disk',radi)));
        mask=conv2(IQ2,se,'same'); mask(mask>0)=1;
    end
    mask=imclose(mask,se); % Tanquem forats petits
end
% imagesc(mask);

% Establir anelles del pealing
pgp=ceil(pg/DX); % Radi pealing en pixels
% masko=mask;
suma=mask;
se = strel('disk',pgp);
se=double(getnhood(se));
areamask=numel(find(mask==1));
while areamask>1,
    
    % erosionar la mask
    % establir si es una cell o un zoom
    
    %     quadrat=zeros(size(Im));
    %     quadrat(1:2,:)=1; quadrat(end-1:end,:)=1;
    %     quadrat(:,1:2)=1; quadrat(:,end-1:end)=1;
    %     sum=quadrat+mask;
    
    %     if find(sum==2), % si es un zoom
    %         if (as(ii)>45)&&(as(ii)<=135),
    %             angle=90;
    %         else
    %             angle=0;
    %         end
    %         se = strel('line',pgp,angle);
    %     else % es tota una cell
    
    %     end
    % mask0=mask;
    
    mask=imerode(mask,se,'same');% imagesc(mask);axis image;
    suma=suma+mask; % imagesc(suma);
    areamask=numel(find(mask==1));
    
    % imagesc(mask+mask2);
end

if zz==1, suma=suma(1:m,round(2/DX)+1:n+round(2/DX)); end
pealing=suma;
Qp=Qi2;
Qp=[Qp zeros(size(Qp,1),1)];

for ii=1:size(Qp,1)
    Qp(ii,3)=pealing(Qi2(ii,2),Qi2(ii,1));
end


end