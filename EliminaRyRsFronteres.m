function [pos2, pos, ind]=EliminaRyRsFronteres(z,Iryr,pos,DX)
% zz   % Zona que elimina els RyRs, entrar en micres. Ex. z=.33;
zz=round(z/DX);

[a,b]=size(Iryr);
% IQ=zeros([a,b]); IQmask=IQ;
% IQmask(1:zz,1:end)=1;IQmask(a-zz:a,1:end)=1;
% IQmask(1:end,1:zz)=1;IQmask(1:end,b-zz:b)=1;
% % imagesc(IQmask);axis image;
% for ii=1:size(pos,1),IQ(pos(ii,1),pos(ii,2))=1;end
% % imagesc(IQ+IQmask);axis image;
% II=IQ+IQmask;
% [a,b]=find(II==2);

ind=[];
for ii=1:size(pos,1)
    if (pos(ii,1)<=zz)||(pos(ii,1)>=a-zz), ind=[ind,ii]; end
    if (pos(ii,2)<=zz)||(pos(ii,2)>=b-zz), ind=[ind,ii]; end
end
 pos2=pos; pos2(ind,:)=[];
 
% % Comprovar:
% IQ=zeros([a,b]);
% IQ2=IQ;
% for ii=1:size(pos,1),IQ(pos(ii,1),pos(ii,2))=1;end
% for ii=1:size(pos2,1),IQ2(pos2(ii,1),pos2(ii,2))=1;end
% % II=IQ+IQ2; % imagesc(II);

end