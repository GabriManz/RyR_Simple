function []=representaIm(Qcom, Inorm, Ch, DX, folder, RF,sr,zoom)

fact=1;
if zoom==0, fact=ceil(.5/DX); end 
cmap=gray(256);
vect=[1:size(cmap,2)];
vect(vect==Ch)=[];
cmap(:,vect)=0; cmap(end,:)=[1 1 1];
if sr==1, cmap(end-1,:)=[1 0 0]; end

ind=find(folder=='\',1,'last');
file=folder(ind+1:end);

for ii=1:length(Qcom)
    filename=['RyRsImage',num2str(ii),'_',file];
    I=zeros(size(Inorm(:,:,1)));
    Qi=Qcom{ii};
    Im=Inorm(:,:,ii); I=imresize(Im,fact,'nearest');I=253*I;
    %(I-min(min(I)))/(max(max(I))-min(min(I)));
   
    for kk=2:length(Qi)
        
        color=255;
        if sr==1,
            if ~strcmp(Qi{kk,9},'Ok'),
                color=254;
            end
            I=textIm(Qi{kk,1}*fact,Qi{kk,2}*fact,num2str(kk-1),I,'textcolor',color,'blending','off');
        else
            if(strcmp(Qi{kk,9},'Ok')),
                I=textIm(Qi{kk,1}*fact,Qi{kk,2}*fact,num2str(kk-1),I,'textcolor',color,'blending','off');
            end
        end
    end
    % imagesc(I);colormap(cmap);axis image;
    
    I=enxufaLlegenda2(I,DX/fact,1);
    
    imwrite(uint8(I),cmap,[folder,'\',RF,'\',filename,'.png']);
    imwrite(uint8(Im*254),cmap,[folder,'\',RF,'\Original_',filename,'.png']);
    
    % Guarda
end

Qi=Qcom{1};
ct=0;
color=255;
I=zeros(size(Inorm(:,:,1)));
Im=Inorm(:,:,ii); I=imresize(Im,fact,'nearest');I=253*I;
for kk=2:length(Qi)
    if strfind(Qi{kk,9},'Ok'),ct=ct+1;
  I=textIm(Qi{kk,1}*fact,Qi{kk,2}*fact,num2str(ct),I,'textcolor',color,'blending','off');
    end
end
imwrite(uint8(I),cmap,[folder,'\',RF,'\Filt_',filename,'.png']);

end

