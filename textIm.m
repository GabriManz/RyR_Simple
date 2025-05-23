function [Im] = textIm(varargin)
%
% outputImage = textIm(x,y,text,image,options)
%
% options can be any of the following couples
%  'fontsize',number (default is 15)
%  'fontname',name (default is 'verdana')
%  'textcolor',[three value vector over unity] (default is [1 1 1])
%  'horizontalalignment', either 'left', 'center' or 'right' (default is 'left')
%  'verticalalignment', either 'top', 'mid' or 'bot' (default is 'mid')
%  'blending', either 'on' or 'off' (default is 'on')
%
% example call:
%     Im2=textIm(100,50,'helloWorld!',originalImage,...
%      'fontsize',8,'fontname','verdana','textcolor',[1 1 0],...
%      'horizontalalignment','center','verticalalignment','bot');

x=varargin{1};
y=varargin{2};
txt=varargin{3};
Im=varargin{4};



colt=[1 1 1];
c=[0.5,0];%centre del text
fs=15;
fn='verdana';
blend=1;
for ii=5:2:nargin
    gaga=lower(varargin{ii});
    gaga=gaga(1:6);
    switch gaga
        case 'fontsi'
            fs=varargin{ii+1};
        case 'textco'
            colt=varargin{ii+1};
        case 'fontna'
            fn=varargin{ii+1}; 
        case 'horizo'
            gaga2=lower(varargin{ii+1}(1:3));
            switch gaga2
                case 'lef'
                    c(2)=0;
                case {'cen','mid'}
                    c(2)=0.5;
                case 'rig'
                    c(2)=1;
            end
        case 'vertic'
           gaga2=lower(varargin{ii+1}(1:3));
            switch gaga2
                case 'top'
                    c(1)=0;
                case {'cen','mid'}
                    c(1)=0.5;
                case 'bot'
                    c(1)=1;
            end
        case 'blendi'
            gaga2=lower(varargin{ii+1}(1:2));
             switch gaga2
                case 'on'
                    blend=1;
                case 'of'
                    blend=0;
             end
        
    end
end

T=1;
if(max(max(max(Im)))<=1),T=1;end
if(max(max(max(Im)))>256),T=max(max(max(Im)));end
if(strcmp(class(Im),'double')&&(T==255)),T=1;end
tenimchars=exist('charset.mat');
if((strcmp(fn,'verdana'))&&(tenimchars==2))
load('charset.mat');


F=[];
for ii=1:length(txt)
    ind=strfind(index,txt(ii));
    if(~isempty(ind))
   F=[F chars{ind(1)}]; 
    end
end


else

figure(424);axes('position',[0,0,1,1]);
plot(0,0,'w+');hold on;plot(1000,1000,'w+');xlim([10 990]);ylim([10 990]);

text(500,500,txt,'fontname',fn,'fontsize',fs);

set(gca,'xtick',[0 1000],'ytick',[0 1000]);
F=getframe(424);F=double(F.cdata);
close(424);

F=F(:,:,1)+F(:,:,2)+F(:,:,3);F=F/3;
F=255-F;F=F/255;
F=F(3:end-3,3:end-3,:);


[a,b]=find(F>0);
F=F(min(a):max(a),min(b):max(b),:);


end

F=imresize(F,fs/size(F,1));
F=(F-min(min(F)))/(max(max(F))-min(min(F)));
if(min(F(1,:))>0)
    F=(F-min(F(1,:)))/(max(max(F))-min(F(1,:)));
    F(F<0)=0;
end

[sv,sh]=size(F);

c=round([sv,sh].*c);


rgy=y-c(1)+1:y+sv-c(1);
rgx=x-c(2)+1:x+sh-c(2);

iny=find((rgy>0)&(rgy<=size(Im,1)));
inx=find((rgx>0)&(rgx<=size(Im,2)));
F=F(iny,inx);
ret=double(Im(rgy(iny),rgx(inx),:));
% figure;imagesc(ret/T)
r2=zeros(size(ret));
placa=ones(size(F));
if(blend==0),F=double(F>.55);end
for ii=1:size(ret,3)
    curr=ret(:,:,ii);
   dest=placa*colt(ii)*T;
   r2(:,:,ii)=curr+(dest-curr).*F;
    
end

% figure;imagesc(F)
% figure;imagesc(r2/T)

if(strcmp(class(Im),'uint8')==1)
r2=uint8(r2);
end
if(strcmp(class(Im),'uint16')==1)
r2=uint16(r2);
end
Im(rgy(iny),rgx(inx),:)=r2;







end