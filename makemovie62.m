function [] = makemovie62(volum,path,DT,tfi,xfi,slo,sho, ch, nom, minmax,cm)

% Fer vídeo d'un volum carregat al workspace:
aux=version;
new=double(str2double(aux(1))>7);

if(nargin<7),sho=1;end % show numeros
if(nargin<6),slo=1;end % slow factor
if(nargin<5),xfi=1;end % space filter
if(nargin<4),tfi=1;end % time filter
if(nargin<3),tag='.';end % id
if(nargin<2),DT=83.333;end 
if(nargin<8),ch=1;end
if(nargin<9),nom='video';end
if(nargin<10),minmax=0;end 
if(nargin<11),cm=gray(256);end

xfilt=ones(xfi,xfi);
xfilt=xfilt/sum(sum(xfilt));
path(path=='/')='\';
if(path(end)=='\'),path=path(1:end-1);end
% nom=path(find(path=='\',2,'last')+1:end);
% nom='video';
sf=tfi;
migsf=sf;
% S2=dir([path '/*' tag '*']);
% if(isempty(S2)),error('No files!');end
% I2=imread([path '/' S2(1).name]);
% [a,b,c]=size(I2);
[a,b,c]=size(volum);

% ch=1;

% if(c>1)
%     aux1=[];
%     for ii=1:c
%         aux1(ii)=sum(sum(I2(:,:,ii)));
%     end
%     [aux,ch]=max(aux1);
% end
% IC=zeros(a,b,sf);
% for ii=1:sf
%     I2=volum(:,:,ii);
%     %I2=I2(:,:,ch);
%     IC(:,:,ii)=I2;
% end
% film=zeros(a-xfi+1,b-xfi+1,c-sf);
% migsf=round(sf/2);
% for ii=sf+1:c
%     for kk=1:sf-1
%         IC(:,:,kk)=IC(:,:,kk+1) ;
%     end
%     I2=volum(:,:,ii);
%     %I2=I2(:,:,ch);
%     IC(:,:,sf)=I2;
%     fr=I2;%sum(IC,3);
%     if(xfi>1)
%         fr=conv2(fr,xfilt,'valid');
%     end
%     film(:,:,ii-sf)=fr;
% end
% q=reshape(film(:,:,1),1,(a-xfi+1)*(b-xfi+1));
% mn=quantile2(q,.001);mx=quantile2(q,.999);
film=volum;
if minmax==1,
    mx=max(max(max(film)));
    mn=min(min(min(film)));
    film=(film-mn)/(mx-mn);
    film(film>1)=1;film(film<0)=0;
end


%  cm=superjet(255);
%  cm(256,:)=[1 1 1];

%  cm=gray(256);
% 
% % l=length(cm);
% for ii=1:3
%     if(ch~=ii)
%         cm(:,ii)=0;
%     end
% end
% %cm(254,:)=[0 0 1];cm(255,:)=[0 1 0];
% cm(end,:)=[1 1 1];
% cm(end-1,:)=[1 0 0];
% %cm(end-5:end-1,3)=1;

if(new==0)
    aviobj = avifile([nom '.avi'],'colormap',cm,'compression','MSVC','fps',round(1000/DT/slo),'quality',100);
else
    aviobj = VideoWriter([nom '.avi'],'Indexed AVI');
    aviobj.FrameRate=round(1000/DT/slo);aviobj.Colormap=cm;
    open(aviobj);
end
for ii=1:size(film,3)
    py=1;
    if(sho==1)
         film(:,:,ii)=textIm(1,py,['Frame ' num2str(ii+migsf-1)],film(:,:,ii),'verticalalignment','top','textcolor',1);
   
        if(slo>1)
            py=py+12;
            film(:,:,ii)=textIm(1,py,['Slow x' num2str(slo)],film(:,:,ii),'verticalalignment','top','textcolor',1);
        end
        if(tfi>1)
            py=py+12;
            film(:,:,ii)=textIm(1,py,['Time mean =' num2str(sf)],film(:,:,ii),'verticalalignment','top','textcolor',1);
        end
        if(xfi>1)
            py=py+12;
            film(:,:,ii)=textIm(1,py,['Space mean =' num2str(xfi)],film(:,:,ii),'verticalalignment','top','textcolor',1);
        end
        
    end
    %imwrite(film(:,:,ii),[folder '\index\' S(jj).name 'f' num2str(ii) '.png']);
    fr=film(:,:,ii);fr=uint8(round(255*fr));
    if(new==0)
        aviobj = addframe(aviobj,fr);
    else
        writeVideo(aviobj,fr);
    end
end
if(new==0)
    aviobj = close(aviobj);
else
    close(aviobj);
end

movefile([nom '.avi'],[path '/' nom '.avi']);

end

