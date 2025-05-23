function [D] = Blob4Diameters(im,gr,he,pl)
% Measure four diameters of centered blob in a square image
% cutting at he normalized height and using gr \mumel wide sections.
% set pl to 1 for a fancy figure

if(nargin<4),pl=0;end
if(nargin<3),he=.5;end
if(nargin<2),gr=1;end
if(nargin<1),error('No image to measure!');end

[a,b]=size(im);
if(a~=b),error('blob image must be square');end
c=round(a/2);
cg=round(gr/2);cg2=gr-cg;
% horitzontal
s1=im(c-cg+1:c+cg2,:);s1=mean(s1,1);
% vertical
s2=im(:,c-cg+1:c+cg2);s2=s2';s2=mean(s2,1);
% diagonal (paral·lel)
imb=zeros(a,b);
for izi=1:a
imb(izi,izi)=1;try imb(izi-cg+1,izi)=1;end;try imb(izi+cg2,izi)=1;end
end
imb=imb(1:a,1:b);imb=imb.*im;
for ii=1:a
    aux=imb(:,ii);aux(aux==0)=[];
   pri(ii)=mean(aux);
   aux=imb(ii,:);aux(aux==0)=[];
   seg(ii)=mean(aux);
end
s3=(pri+seg)/2;s3=s3(gr:end-gr+1);s3(isnan(s3))=0;

% diagonal (meridiana)
imb=zeros(a,b);
for izi=1:a
imb(a-izi+1,izi)=1;try imb(a-izi+1-cg+1,izi)=1;end;try imb(a-izi+1+cg2,izi)=1;end
end
imb=imb(1:a,1:b);imb=imb.*im;
for ii=1:a
    aux=imb(:,ii);aux(aux==0)=[];
   pri(ii)=mean(aux);
   aux=imb(a-ii+1,:);aux(aux==0)=[];
   seg(ii)=mean(aux);
end
s4=(pri+seg)/2;s4=s4(gr:end-gr+1);s4(isnan(s4))=0;

%figure;plot(s1);hold on;plot(s2,'c');plot(s3,'r');plot(s4,'m');
[ft1,rsquare1] = blobfit(s1');[D1] = tallaDiametre(ft1,he);
[ft2,rsquare2] = blobfit(s2');[D2] = tallaDiametre(ft2,he);
[ft3,rsquare2] = blobfit(s3');[D3] = tallaDiametre(ft3,he);
[ft4,rsquare2] = blobfit(s4');[D4] = tallaDiametre(ft4,he);
D=mean([D1,D2,D3*sqrt(2),D4*sqrt(2)]);

if(pl==1),FWHM(im,gr,he,DX);end

end


function [ft,rsquare] = blobfit(s)
s=(s-min(s))/(max(s)-min(s));
len=length(s);
c=len/2;
% exponencial simple

 model1 = fittype('exp(-((x-m)^2)/(2*(s^2)))');
        
            options1 = fitoptions(model1);
            options1.Display = 'off';
            options1.Robust = 'bisquare';
            %%%%%%%%%%%%%%%%%%%    m s     
            options1.startpoint = [c 1];
            options1.Lower = [0 0];
            options1.Upper = [length(s) length(s)/2];
            [fresult1,gof2,out2]  = fit([1:len]',s,model1,options1);
            rsquare=gof2.rsquare;

      
ft=exp(-(([1:len]-fresult1.m).*([1:len]-fresult1.m))/(2*(fresult1.s^2)));



end

function [D] = tallaDiametre(ft,he)

fti=interp(ft,10);
x=[10:length(fti)+9]/10;
gg=fti-he;
gg=gg(1:end-1).*gg(2:end);
o1=find(gg<0,1,'first');o2=find(gg<0,1,'last')+1;
o1=x(o1);o2=x(o2);
D=abs(o1-o2);


end
function [D,o1,o2] = tallaDiametre2(ft,he)

fti=interp(ft,10);
x=[10:length(fti)+9]/10;
gg=fti-he;
gg=gg(1:end-1).*gg(2:end);
o1=find(gg<0,1,'first');o2=find(gg<0,1,'last')+1;
oo1=x(o1);oo2=x(o2);
D=abs(oo1-oo2);


end

function [Dmean] = FWHM(fat,gr,he,DX)
close all;

%%%%%%%%%%%%%%%%%
% general stuff %
%%%%%%%%%%%%%%%%%
FS=28;% fontsize
greenCM=gray(255);greenCM(:,2)=0;greenCM(:,3)=0;cols=superjet(255);
cmap=greenCM;% image colormap
%cols=cols([175,194,233,243],:);% gay 
%cols=cols([175,194,228,233],:);% whine
%cols=cols([175,194,55,73],:);% blues
cols=cols([96,156,55,73],:);% blues
fac=1;% color factor for fit related to diameter
fac2=.3;% color factor for raw related to diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D panel & zenital view %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
NF1=1;% normalizing factor
x_l='X(\mum)';y_l='Y(\mum)';% axis labels
st1='-';st2='*';% linestyles for fit and diameter line
ms1=6;%marker size in diameter line in plot 3D
ms2=8;%marker size in diameter line in zenital view
wi1=2;wi2=2;% line widths for fit and diameter line
%%%%%%%%%%%%%%
% 4subpanels %
%%%%%%%%%%%%%%
boxonoff='off';% show box in axes
wi3=4;wi4=1;% linewidths for fit and plot
wi5=3;% linewidths for diameter
st3='-';st4='o';st5='-+';% linestyles for fit and plot
pt='lin';% plot type (lin or bar)
ox=0.3333;% x position offset
tv=.45;th=.16;% vertical and horizontal sizes
mv=.045;mh=.02;% vertical and horizontal margins
NF2=0;%Normalizing factor in 
ylyl=[min(min(fat)),max(max(fat))];% ylimits when NF2~=1
ylylalt=[-0.05 1.05];% ylimits when NF2=1
fwtextX=0.5;%ylyl(1)+diff(ylyl)*.9;% text in X position
FL='F';% ylabel when NF2~=1
FLalt='F';% ylabel when NF2=1;
% 
% load([RF 'zMetaData.mat']);
% load([RF 'zData0.mat']);
% load([RF 'zData2.mat']);


% DX=1;

  %try  

[a,b]=size(fat);
if(a~=b),error('blob image must be square');end
c=round(a/2);
cg=round(gr/2);cg2=gr-cg;
% horitzontal
s1=fat(c-cg+1:c+cg2,:);s1=mean(s1,1);
% vertical
s2=fat(:,c-cg+1:c+cg2);s2=s2';s2=mean(s2,1);
% diagonal (paral·lel)
imb=zeros(a,b);
for izi=1:a
imb(izi,izi)=1;try imb(izi-cg+1,izi)=1;end;try imb(izi+cg2,izi)=1;end
end
imb=imb(1:a,1:b);imb=imb.*fat;
for izi=1:a
    aux=imb(:,izi);aux(aux==0)=[];
   pri(izi)=mean(aux);
   aux=imb(izi,:);aux(aux==0)=[];
   seg(izi)=mean(aux);
end
s3=(pri+seg)/2;s3=s3(gr:end-gr+1);s3(isnan(s3))=0;

% diagonal (meridiana)
imb=zeros(a,b);
for izi=1:a
imb(a-izi+1,izi)=1;try imb(a-izi+1-cg+1,izi)=1;end;try imb(a-izi+1+cg2,izi)=1;end
end
imb=imb(1:a,1:b);imb=imb.*fat;
for izi=1:a
    aux=imb(:,izi);aux(aux==0)=[];
   pri(izi)=mean(aux);
   aux=imb(a-izi+1,:);aux(aux==0)=[];
   seg(izi)=mean(aux);
end
s4=(pri+seg)/2;s4=s4(gr:end-gr+1);s4(isnan(s4))=0;

%figure;plot(s1);hold on;plot(s2,'c');plot(s3,'r');plot(s4,'m');
[ft1,rsquare1] = blobfit(s1');[D1,o1,f1] = tallaDiametre2(ft1,he);
[ft2,rsquare2] = blobfit(s2');[D2,o2,f2] = tallaDiametre2(ft2,he);
[ft3,rsquare2] = blobfit(s3');[D3,o3,f3] = tallaDiametre2(ft3,he);
[ft4,rsquare2] = blobfit(s4');[D4,o4,f4] = tallaDiametre2(ft4,he);
Dmean=mean([D1,D2,D3*sqrt(2),D4*sqrt(2)])*DX;
D3=D3*sqrt(2);
D4=D4*sqrt(2);



numeret=100*gr;

figure(numeret); 
set(gcf,'color','w','position',[350   784   1800   504]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% panell amb imatge i gaussianes 3D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes('position',[.666 0 .333 1.2],'fontsize',FS);
hold on;
view(37,26);
view(38,29);
surface('XData',[0 a-1;0 a-1]*DX,'YData',[0 0;b-1 b-1]*DX,...
        'ZData',[0 0; 0 0],'CData',flipdim(fat,1),...
        'FaceColor','texturemap','EdgeColor','none');
    colormap(cmap);
    zlim([0,1]); pause(.1);  
xtl=get(gca,'xticklabel'); 
xl=xlabel(x_l);set(xl,'position',[(size(fat,2)*0.5*DX),(size(fat,2)*-0.1*DX),0]);
yl=ylabel(y_l);set(yl,'position',[(size(fat,2)*1.1*DX),(size(fat,2)*0.5*DX),0]); 


if(NF1~=1),NF1=max(s1);end
fti=interp(ft1,10);fti=fti-min(fti); x=[10:length(fti)+9]/10;x=x-1;
h1=plot3(x(1:310)*DX,(DX*(a-1)/2)*ones(1,310),NF1*fti(1:310));
set(h1,'linewidth',wi1,'linestyle',st1,'color',cols(1,:)*fac);
d1=plot3([x(o1) x(f1)]*DX,[(DX*(a-1)/2) (DX*(a-1)/2)],NF1*[fti(o1) fti(f1)],st2);
set(d1,'linewidth',wi2,'color',cols(1,:),'markersize',ms1);
if(NF1~=1),NF1=max(s2);end
fti=interp(ft2,10);fti=fti-min(fti);
y=[10:length(fti)+9]/10;y=y-1;y=y(1:310);y=y(end:-1:1);fti=fti(1:310);
h2=plot3((DX*(b-1)/2)*ones(1,310),y*DX,NF1*fti);
set(h2,'linewidth',wi1,'linestyle',st1,'color',cols(2,:)*fac);
d2=plot3([(DX*(a-1)/2) (DX*(a-1)/2)],[y(o2) y(f2)]*DX,NF1*[fti(o2) fti(f2)],st2);set(d2,'linewidth',wi2,'color',cols(2,:),'markersize',ms1);
if(NF1~=1),NF1=max(s3);end
fti=interp(ft3,10); fti=fti-min(fti);
x=[1:length(ft3)];x=x+gr-1;x=x(1:length(ft3));y=x(end:-1:1);
x=interp(x,10);x=x-1;
y=interp(y,10);y=y-1;
h3=plot3(x*DX,y*DX,NF1*fti);
set(h3,'linewidth',wi1,'linestyle',st1,'color',cols(3,:)*fac);
d3=plot3([x(o3) x(f3)]*DX,[y(o3) y(f3)]*DX,NF1*[fti(o3) fti(f3)],st2);
set(d3,'linewidth',wi2,'color',cols(3,:),'markersize',ms1);
if(NF1~=1),NF1=max(s4);end
fti=interp(ft4,10);fti=fti-min(fti);
x=[1:length(ft4)];x=x+gr-1;x=x(1:length(ft3));y=x;
x=interp(x,10);x=x-1;
y=interp(y,10);y=y-1;
h4=plot3(x*DX,y*DX,NF1*fti);set(h4,'linewidth',wi1,'linestyle',st1,'color',cols(4,:)*fac);
d4=plot3([x(o4) x(f4)]*DX,[y(o4) y(f4)]*DX,NF1*[fti(o4) fti(f4)],st2);set(d4,'linewidth',wi2,'color',cols(4,:),'markersize',ms1);


% if(NF1~=1),NF1=max(s1);end
% fti=interp(ft1,10);x=[10:length(fti)+9]/10;x=x-1;
% h1=plot3(x(:)*DX,(DX*(a-1)/2)*ones(1,310),NF1*fti(1:310));
% set(h1,'linewidth',wi1,'linestyle',st1,'color',cols(1,:)*fac);
% d1=plot3([x(o1) x(f1)]*DX,[(DX*(a-1)/2) (DX*(a-1)/2)],NF1*[fti(o1) fti(f1)],st2);
% set(d1,'linewidth',wi2,'color',cols(1,:),'markersize',ms1);
% if(NF1~=1),NF1=max(s2);end
% fti=interp(ft2,10);y=[10:length(fti)+9]/10;
% y=y-1;y=y(1:310);y=y(end:-1:1);fti=fti(1:310);
% h2=plot3((DX*(b-1)/2)*ones(1,310),y*DX,NF1*fti);
% set(h2,'linewidth',wi1,'linestyle',st1,'color',cols(2,:)*fac);
% d2=plot3([(DX*(a-1)/2) (DX*(a-1)/2)],[y(o2) y(f2)]*DX,NF1*[fti(o2) fti(f2)],st2);
% set(d2,'linewidth',wi2,'color',cols(2,:),'markersize',ms1);
% if(NF1~=1),NF1=max(s3);end
% fti=interp(ft3,10);
% x=[1:length(ft3)];x=x+gr-1;x=x(1:length(ft3));y=x(end:-1:1);
% x=interp(x,10);x=x-1;
% y=interp(y,10);y=y-1;
% h3=plot3(x*DX,y*DX,NF1*fti);set(h3,'linewidth',wi1,'linestyle',st1,'color',cols(3,:)*fac);
% d3=plot3([x(o3) x(f3)]*DX,[y(o3) y(f3)]*DX,NF1*[fti(o3) fti(f3)],st2);set(d3,'linewidth',wi2,'color',cols(3,:),'markersize',ms1);
% if(NF1~=1),NF1=max(s4);end
% fti=interp(ft4,10);
% x=[1:length(ft4)];x=x+gr-1;x=x(1:length(ft3));y=x;
% x=interp(x,10);x=x-1;
% y=interp(y,10);y=y-1;
% h4=plot3(x*DX,y*DX,NF1*fti);set(h4,'linewidth',wi1,'linestyle',st1,'color',cols(4,:)*fac);
% d4=plot3([x(o4) x(f4)]*DX,[y(o4) y(f4)]*DX,NF1*[fti(o4) fti(f4)],st2);set(d4,'linewidth',wi2,'color',cols(4,:),'markersize',ms1);


set(gca,'ztick',[],'ZColor','w');axis square;
zlim([0,NF1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% panell amb imatge i diametres 2D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=axes('position',[0.01 .05 .33 .93],'fontsize',FS);
hold on;
view(0,90);
surface('XData',[0 a-1;0 a-1]*DX,'YData',[0 0;b-1 b-1]*DX,...
        'ZData',[0 0; 0 0],'CData',flipdim(fat,1),...
        'FaceColor','texturemap','EdgeColor','none');
    colormap(greenCM)
    zlim([0,1]);
xtl=get(gca,'xticklabel'); 
xl=xlabel(x_l);set(xl,'position',[str2num(char(xtl(round(size(xtl,1)/2),:))),-size(fat,1)/60],'backgroundcolor','w');
ylabel(y_l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% panell amb 4 subpanells i gaussianes 1D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B1=axes('position',[ox 2*mv+tv th tv],'fontsize',FS);
B2=axes('position',[ox+th+mh 2*mv+tv th tv],'fontsize',FS);
B3=axes('position',[ox mv th tv],'fontsize',FS);
B4=axes('position',[ox+th+mh mv th tv],'fontsize',FS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pinta plots, fits i linies %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wi1=3;
axes(A);
if(NF2==1),s1=n1D(s1);ylyl=ylylalt;FL=FLalt;end
fti=interp(ft1,10);x=[10:length(fti)+9]/10;x=x-1;
h1=plot3(x(1:310)*DX,(DX*(a-1)/2)*ones(1,310),fti(1:310));set(h1,'linewidth',wi1,'linestyle',st1,'color',cols(1,:)*fac);
d1=plot3([x(o1) x(f1)]*DX,[(DX*(a-1)/2) (DX*(a-1)/2)],10*[fti(o1) fti(f1)],st2);set(d1,'linewidth',wi2,'color',cols(1,:),'markersize',ms2);
axes(B1);
fti2=((max(s1)-min(s1))*fti(1:310))+min(s1);
if(strcmp(pt,'lin')==1)
p1=plot([0:length(s1)-1]*DX,s1,st4);set(p1,'linewidth',wi4,'color',cols(1,:)*fac2);hold on;
end
if(strcmp(pt,'bar')==1)
p1=bar([0:length(s1)-1]*DX,s1);set(p1,'linewidth',wi4,'facecolor',cols(1,:)*fac2);hold on;
end
if(strcmp(pt,'ste')==1)
p1=stem([0:length(s1)-1]*DX,s1);set(p1,'linewidth',wi4,'color',cols(1,:)*fac2);hold on;
end
q1=plot([x(1:310)]*DX,fti2);set(q1,'linewidth',wi3,'linestyle',st3,'color',cols(1,:)*fac);
d1=plot(([x(o1) x(f1)])*DX,[fti2(o1) fti2(f1)],st5);set(d1,'linewidth',wi5,'color',cols(1,:));
axis tight;ylim(ylyl);ylabel(FL);xlxl=xlim();h=text(xlxl(1)+fwtextX*diff(xlxl),ylyl(1)+diff(ylyl)*.9,['FWHM_1 = ' num2str(0.1*round(D1*DX*10)) '\mum']);set(h,'fontsize',FS,'horizontalalignment','center');
set(gca,'box',boxonoff);
axes(A);
if(NF2==1),s2=n1D(s2);end
fti=interp(ft2,10);y=[10:length(fti)+9]/10;y=y-1;y=y(1:310);y=y(end:-1:1);fti=fti(1:310);
h2=plot3((DX*(b-1)/2)*ones(1,310),y*DX,fti);set(h2,'linewidth',wi1,'linestyle',st1,'color',cols(2,:)*fac);
d2=plot3([(DX*(a-1)/2) (DX*(a-1)/2)],[y(o2) y(f2)]*DX,10*[fti(o2) fti(f2)],st2);set(d2,'linewidth',wi2,'color',cols(2,:),'markersize',ms2);
axes(B2);
fti2=((max(s2)-min(s2))*fti(1:310))+min(s2);
if(strcmp(pt,'lin')==1)
p2=plot([0:length(s2)-1]*DX,s2,st4);set(p2,'linewidth',wi4,'color',cols(2,:)*fac2);hold on;
end
if(strcmp(pt,'bar')==1)
p2=bar([0:length(s2)-1]*DX,s2);set(p2,'linewidth',wi4,'facecolor',cols(2,:)*fac2);hold on;    
end
q2=plot([x(1:310)]*DX,fti2);set(q2,'linewidth',wi3,'linestyle',st3,'color',cols(2,:)*fac);
d2=plot(([x(o2) x(f2)])*DX,[fti2(o2) fti2(f2)],st5);set(d2,'linewidth',wi5,'color',cols(2,:));
axis tight;ylim(ylyl);xlxl=xlim();h=text(xlxl(1)+fwtextX*diff(xlxl),ylyl(1)+diff(ylyl)*.9,['FWHM_2 = ' num2str(0.1*round(D2*DX*10)) '\mum']);set(h,'fontsize',FS,'horizontalalignment','center');
set(gca,'box',boxonoff);
axes(A);
if(NF2==1),s3=n1D(s3);end
fti=interp(ft3,10);
x=[1:length(ft3)];x=x+gr-1.5;x=x(1:length(ft3));y=x(end:-1:1);
x=interp(x,10);x=x-1;
y=interp(y,10);%y=y+1;
h3=plot3(x*DX,y*DX,fti);set(h3,'linewidth',wi1,'linestyle',st1,'color',cols(3,:)*fac);
d3=plot3([x(o3) x(f3)]*DX,[y(o3) y(f3)]*DX,10*[fti(o3) fti(f3)],st2);set(d3,'linewidth',wi2,'color',cols(3,:),'markersize',ms2);
axes(B3);
fti2=((max(s3)-min(s3))*fti)+min(s3);
if(strcmp(pt,'lin')==1)
p3=plot((gr-1+[1:length(s3)])*sqrt(2)*DX,s3,st4);set(p3,'linewidth',wi4,'color',cols(3,:)*fac2);hold on;
end
if(strcmp(pt,'bar')==1)
p3=bar((gr-1+[1:length(s3)])*sqrt(2)*DX,s3);set(p3,'linewidth',wi4,'facecolor',cols(3,:)*fac2);hold on;    
end
q3=plot([1+x]*DX*sqrt(2),fti2);set(q3,'linewidth',wi3,'linestyle',st3,'color',cols(3,:)*fac);
d3=plot((1+[x(o3) x(f3)])*DX*sqrt(2),[fti2(o3) fti2(f3)],st5);set(d3,'linewidth',wi5,'color',cols(3,:));
axis tight;ylim(ylyl);ylabel(FL);
xtl=get(gca,'xticklabel'); xl=xlabel('d(\mum)');set(xl,'position',[.9*length(s3)*sqrt(2),-diff(ylyl)*.095],'horizontalalignment','center','verticalalignment','middle');%set(xl,'position',[(str2num(xtl(end,:))+str2num(xtl(end-1,:)))/2,-diff(ylyl)*.095],'horizontalalignment','center','verticalalignment','middle');
xlxl=xlim();h=text(xlxl(1)+fwtextX*diff(xlxl),ylyl(1)+diff(ylyl)*.9,['FWHM_3 = ' num2str(0.1*round(D3*DX*10)) '\mum']);set(h,'fontsize',FS,'horizontalalignment','center');
set(gca,'box',boxonoff);
axes(A);
if(NF2==1),s4=n1D(s4);end
fti=interp(ft4,10);
x=[1:length(ft4)];x=x+gr-1.5;x=x(1:length(ft3));y=x;
x=interp(x,10);x=x-1;
y=interp(y,10);y=y-1;
h4=plot3(x*DX,y*DX,fti);set(h4,'linewidth',wi1,'linestyle',st1,'color',cols(4,:)*fac);
d4=plot3([x(o4) x(f4)]*DX,[y(o4) y(f4)]*DX,10*[fti(o4) fti(f4)],st2);set(d4,'linewidth',wi2,'color',cols(4,:),'markersize',ms2);
axes(B4);
fti2=((max(s4)-min(s4))*fti)+min(s4);
if(strcmp(pt,'lin')==1)
p4=plot((gr-1+[1:length(s4)])*sqrt(2)*DX,s4,st4);set(p4,'linewidth',wi4,'color',cols(4,:)*fac2);hold on;
end
if(strcmp(pt,'bar')==1)
p4=bar((gr-1+[1:length(s4)])*sqrt(2)*DX,s4);set(p4,'linewidth',wi4,'facecolor',cols(4,:)*fac2);hold on;    
end
q4=plot([1+x]*DX*sqrt(2),fti2);set(q4,'linewidth',wi3,'linestyle',st3,'color',cols(4,:)*fac);
d4=plot((1+[x(o4) x(f4)])*DX*sqrt(2),[fti2(o4) fti2(f4)],st5);set(d4,'linewidth',wi5,'color',cols(4,:));
axis tight;ylim(ylyl);%ylabel('{\Delta}F/F_0');
xtl=get(gca,'xticklabel'); xl=xlabel('d(\mum)');set(xl,'position',[.9*length(s4)*sqrt(2),-diff(ylyl)*.095],'horizontalalignment','center','verticalalignment','middle');
xlxl=xlim();h=text(xlxl(1)+fwtextX*diff(xlxl),ylyl(1)+diff(ylyl)*.9,['FWHM_4 = ' num2str(0.1*round(D4*DX*10)) '\mum']);set(h,'fontsize',FS,'horizontalalignment','center');
set(gca,'box',boxonoff);
axes(A);
 set(gca,'ztick',[]);axis square;
% LE=legend([d1,d2,d3,d4],['FWHM_1 = ' num2str(0.1*round(D1*DX*10)) '\mum'],['FWHM_2 = ' num2str(0.1*round(D2*DX*10)) '\mum'],['FWHM_3 = ' num2str(0.1*round(D3*DX*10)) '\mum'],['FWHM_4 = ' num2str(0.1*round(D4*DX*10)) '\mum']);
% LP=get(LE,'position');LP(1)=LP(1)+.53;set(LE,'position',LP);
% 
% LP(4)=(LP(4)/4);
% LP(2)=LP(2)-LP(4)-.01;
annotation('textbox',[.83 .7 .2 .1],'fontsize',FS,'horizontalalignment','center','edgecolor','w','string',['<FWHM> = ' num2str(0.01*round(mean([D1,D2,D3,D4])*DX*100)) '\mum']);
% LP=get(LE,'position');LP(1)=LP(1)-.18;
% annotation('textbox',[LP],'fontsize',FS,'horizontalalignment','left','edgecolor','w','string',['Section thickness:' num2str(gr) '\mum']);


saveWysiwyg(gcf,(['D:\Dades Lab\RyR-GFP_Activation_Sparks\FIGURAS\190411\SPK_FWHM_28.png']));

%[ft1,rsquare1] = blobfit(s1');[D1,o1,f1] = tallaDiametre(ft1,he);
  % end


end