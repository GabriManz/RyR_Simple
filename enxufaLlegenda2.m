function [Im] = enxufaLlegenda2(Im,res,col,pos,sty)
% Incrusta barra d'escala tenint en compte que
% cada pixel fa res micres. L'escala es pinta amb color col i estil sty 
% (0 per linia simple, 1 per linia amb extrems, 2 per barra).
% pos pot ser un de 'NE','NW','SE','SW' (default es 'SW')





if(nargin<5),sty=1;end
if(nargin>=4),
    if(strcmp(pos,'SE')==1),pos=1;end
    if(strcmp(pos,'SW')==1),pos=2;end
    if(strcmp(pos,'NE')==1),pos=3;end
    if(strcmp(pos,'NW')==1),pos=4;end
    if(ischar(pos)),pos=1;end
end
if(nargin<4),pos=1;end
if(nargin<3),error('falten arguments');end

T=255;
if(max(max(max(Im)))<=1),T=1;end
if(max(max(max(Im)))>256),T=max(max(max(Im)));end

    if(strcmp(class(Im),'uint8')==1),cla='uint8';end
    if(strcmp(class(Im),'uint16')==1),cla='uint16';T=2^16;end
    if(strcmp(class(Im),'double')==1),cla='double';end

Im=double(Im);
[a,b,c]=size(Im);



tamany=[a,b].*[res,res];%en micres

marge=round([b/10 b/5]);%enpix

accepted=[1,2,5,10,20,50,100,200,500,1000,5000,10000,50000,100000]/res;
ind=find((accepted<=marge(2))&(accepted>=marge(1)),1,'first');

longitud=round(accepted(ind));
%F=ceil(a/20);
G=ceil(a/200);
if(pos>2)% a dalt
barraY=3*G:4*G;
palsY=G:6*G;
Factoret=-1;
else  % a baix
barraY=a-4*G:a-3*G;
palsY=a-6*G:a-G;
Factoret=1;
end  
if(mod(pos,2)==1) % dreta
barraX=b-longitud-2*G:b-2*G;
palsX=[b-longitud-2*G:b-longitud-G,b-3*G:b-2*G];
else     % esquerra
barraX=2*G:2*G+longitud;
palsX=[2*G:3*G,longitud+G:longitud+2*G];
end




for ii=1:c
if(sty<=1)
Im(barraY,barraX,ii)=col(ii)*T;
end
if(sty>=1)
Im(palsY,palsX,ii)=col(ii)*T; 
end
if(sty==2)
 Im(palsY(1):palsY(1)+numel(barraY)-1,barraX,ii)=col(ii)*T;
 Im(palsY(end)-numel(barraY)+1:palsY(end),barraX,ii)=col(ii)*T;  
 gru=round(numel(palsX)/2);
 inc=longitud/round(accepted(ind)*res);
 for jj=1:round(accepted(ind)*res)-1
     Im(palsY,palsX(1)+inc*jj:palsX(1)+inc*jj+gru,ii)=col(ii)*T; 
     
 end
end
end
centre=round([mean(barraY) mean(barraX)]);
centre(1)=centre(1)-G*Factoret;

numeret=num2str(round(accepted(ind)*res));
numeret=num2str(numeret);
%Z=text2im([num2str(numeret) '?m']);

load('charset.mat');

Z=[];
for ii=1:length(numeret)
    ind=strfind(index,numeret(ii));
    if(~isempty(ind))
   Z=[Z chars{ind(1)}]; 
    end
end
  Z=[Z chars{145} chars{53}];  
  Z=imresize(Z,ceil(b/10/longitud));
  mn=median([Z(1,:),Z(end,:)]);
  Z=(Z-mn)/(max(max(Z))-mn);
  Z(Z<0)=0;
  if(c==1),Z=round(Z);end
if(pos<=2)  
   centre(1)=centre(1)-size(Z,1)*Factoret;
end
   centre(2)=centre(2)-ceil(size(Z,2)/2);
   ret=Im(centre(1):centre(1)+size(Z,1)-1,centre(2):centre(2)+size(Z,2)-1,:);
   placa=ones(size(Z));
    for ii=1:c
    curr=ret(:,:,ii);
    dest=placa*col(ii)*T;
Im(centre(1):centre(1)+size(Z,1)-1,centre(2):centre(2)+size(Z,2)-1,ii)=curr+(dest-curr).*Z;
    end
   
 eval(['Im=' cla '(Im);']);  
   
%figure;imagesc(1);colormap([0 0 0]);
%text(1,1,[ '1 2 3 4 5 6 7 8 9 0 \mum'],'color',[1 1 1],'fontsize',14,'HorizontalAlignment','center');









end
% 
% % % 
% alph='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
% for ii=1:length(alph)
% %for ii=0:11
%     figure(23);imagesc(1);colormap([0 0 0]);
% % if(ii<=9)
% %     locualo=num2str(ii);
% % end
% % if(ii==10),locualo='\mu';end
% % if(ii==11),locualo='m';end
% locualo=alph(ii);
% text(1,1,locualo,'color',[1 1 1],'fontsize',14,'HorizontalAlignment','center');
%     [aa,bb]=getframe(23);close(23); 
%     aa=rgb2gray(aa);
%     %if(ii==10), aa=aa(239:256,343:352);else
%     [ka,kb]=find(aa>=247);
%     aa=aa(min(ka):max(ka),min(kb):max(kb));
%     %end
%    % if(ii==11), aa=aa(1:end-1,:);aa=[zeros(1,size(aa,2));aa];end
%     aa=aa/max(max(aa));
%     [ga,gb]=size(aa);
%     eny=18-ga;
%     enx=10-gb;
%     for jj=1:eny
%        if(mod(jj,2)==0)
%           aa=[aa;zeros(1,size(aa,2))]; 
%        else
%           aa=[zeros(1,size(aa,2));aa]; 
%        end
%     end
% %     for jj=1:enx
% %        if(mod(jj,2)==0)
%            aa=[aa,zeros(size(aa,1),1)]; 
% %        else
% %           aa=[zeros(size(aa,1),1),aa]; 
% %        end
% %     end
% %imwrite(aa,gray(numel(unique(aa))),['charSet\' num2str(ii) '.png']);
% imwrite(aa,gray(numel(unique(aa))),['charSet\' alph(ii) num2str(ceil(ii/26)) '.png']);
% end

    
    