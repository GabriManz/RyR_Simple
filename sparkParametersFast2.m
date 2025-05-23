function [sortida r2 I bline centre] = sparkParametersFast2(imsize,x,T,centre,plt,titol,cutout,normVal,cols,DT,nnr)
%if(nargin<4),plt=0;else plt=1;end
if(nargin<9),txtfunc=1;else txtfunc=2;end
I=0;
imsize(2)=imsize(2)-1;
% sortida=[amp,tau,ror,t2p,dhm,AMP];
xfp=x;Tfp=T;
if(cutout<length(T))
    x=x(1:cutout);T=T(1:cutout);
end
gru=2;% gruix linies
FZ=14;FS=16;if(txtfunc~=1),FS=FS-1;end
bline=1;
sortida=[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
r2=0;
try
    % AMP
    win=ceil(15/DT);
    vv=max([1 centre(1)-win]):min([nnr centre(1)+win*2]);
    [amp ind]=max(T(vv));% (F-Fo)/Fo
    centre(1)=vv(ind);
    sortida(6)=amp;
    
    %     wps=round(25/DT);
    %     wbl=round(200/DT);
    %     wpost=round(300/DT);
    
    aux1=centre(1);
    t=centre(1);
    sig=T;
    %     in=max([1 t-wbl]);v=[in:t];
    v=[1:t];
    sw=sig(v);
    %     if length(in:t)>6,
    %         quarticMA = sgolayfilt(sig, 8, 19);
    %         swf=quarticMA(in:t); % https://es.wikipedia.org/wiki/Filtro_de_Savitzky%E2%80%93Golay
    %         pend=zeros(1,length(sw)-1); pendf=pend;
    %         for jj=2:length(sw)-1
    %             pend(jj)=(sw(jj)-sw(1))/jj;
    %             pendf(jj)=(swf(jj)-swf(1))/jj;
    %         end
    %         [MM,mm]=findPeaks6(pendf);
    %         if ~isempty(mm),
    %             if length(mm)>1,
    %             if pendf(mm(end))>pendf(mm(end-1)),inblsw=mm(end-1); else, inblsw=mm(end); end
    %             inbl=v(inblsw)+1;
    %             else
    %                 inblsw=mm; inbl=v(inblsw)+1;
    %             end
    %         else
    %             [~,inbl]=min(sw);
    %         end
    %
    %     else
    %         [~,inbl]=min(sw);
    %     end
    swd=diff(diff(sw));
    [MM,mm]=findPeaks6(swd);
    [~,cc]=max(swd(MM));
    inbl=MM(cc)+2;
    aux=inbl;
    
    % %%%%%%%%%%%%%%
    T1=T(aux:aux1);
    x1=x(aux:aux1);
    %figure;plot(x,T);hold on;plot(x1,T1,'r');
    
    % fit lineal a la pujada
    if(size(x1,1)<size(x1,2)),x1=x1';end
    if(size(T1,1)<size(T1,2)),T1=T1';end
    [fresult2,~]=polyfit(x1,T1,1);
    T1=fresult2(1)*x1+fresult2(2);
    fitted1=fresult2(1)*x+fresult2(2);
    %figure;plot(x,fitted1);hold on;plot(x,T);plot(x,fitted2);
    
    % AMP & BL
    bline=quantile(T(1:aux),0.5);
    amp=amp-bline;
    sortida(1)=amp;
    
    % t2p
    %aux=find(T(1:centre(1))<1,1,'last');
    xb=(bline-fresult2(2))/fresult2(1);
    xp=x(aux+round((centre(1)-aux)/2));
    t2p=x(centre(1))-xb;
    sortida(4)=t2p;
    
    
    % RoR
    ror=fresult2(1);
    sortida(3)=ror;
    
    % Tau
    caiguda=T(centre(1):end);
    TT=caiguda;
    xx=x(centre(1):end);
    model1 = fittype('a*exp(-b*(x-e))+f');
    options1 = fitoptions(model1);
    options1.Display = 'off';
    options1.Robust = 'bisquare';
    %%%%%%%%%%%%%%%%%%%%%%      a       b    e        f
    options1.startpoint = [max(abs(TT)) 0  min(xx) min(TT)];
    options1.Lower = [0 0 -10*max(xx) min(TT)];
    options1.Upper = [100*max(abs(TT)) 1 100*max(xx) max((TT))];
    [fresult1,gof2,out2]  = fit(xx',TT',model1,options1);
    rsquare1=gof2.rsquare;
    r2=rsquare1;
    if r2<0.3,  % Dins d'una wave? % By Carmeta!!!!
        win=ceil(100/DT); % 100 ms
        TT=TT(1:win); xx=xx(1:win); % Escursa
        options1.startpoint = [max(abs(TT)) 0  min(xx) min(TT)];
        options1.Lower = [0 0 -10*max(xx) min(TT)];
        options1.Upper = [100*max(abs(TT)) 1 100*max(xx) max((TT))];
        [fresult1,gof2,out2]  = fit(xx',TT',model1,options1);
        rsquare1=gof2.rsquare;
        r2=rsquare1;
    end % Dins d'una wave? % By Carmeta!!!!
    if r2<0.1,  % Dins d'una wave? % By Carmeta!!!!
        win=ceil(60/DT); % 100 ms
        TT=TT(1:win); xx=xx(1:win); % Escursa
        options1.startpoint = [max(abs(TT)) 0  min(xx) min(TT)];
        options1.Lower = [0 0 -10*max(xx) min(TT)];
        options1.Upper = [100*max(abs(TT)) 1 100*max(xx) max((TT))];
        [fresult1,gof2,out2]  = fit(xx',TT',model1,options1);
        rsquare1=gof2.rsquare;
        r2=rsquare1;
    end % Dins d'una wave? % By Carmeta!!!!
    %     figure;plot(x,T,'-b'); hold on;plot(fresult1,'r');
    % figure;plot(xx,TT,'-b'); hold on;plot(fresult1,'r');
    fitted2=fresult1.a*exp(-fresult1.b*(x-fresult1.e))+fresult1.f;
    T2=fresult1.a*exp(-fresult1.b*(xx-fresult1.e))+fresult1.f;
    % figure;plot(x,T,'-b'); hold on;plot(fresult1,'r');
    tau=1/fresult1.b;
    sortida(2)=tau;
    %sortida=[amp,tau,ror,t2p,dhm];
    
    % FDHM
    % AMP50=((amp-bline)/2)+bline; % CANVI 190708!
    AMP50=((sortida(6)-bline)/2)+bline;
    AMP25=amp-((sortida(6)-bline)/4)+bline;
    AMP75=amp-((sortida(6)-bline)*(3/4))+bline;
    AMP90=amp-((sortida(6)-bline)*(90/100))+bline;
    % fitted1 --> RoR
    % fitted2 --> tau
    %     x3=linspace(x1(1),x1(2),length(x1)*4);
    %     xx2=linspace(xx(1),xx(2),length(xx)*4);
    %     T1i=fresult2(1)*x3+fresult2(2);
    %     T2i=fresult1.a*exp(-fresult1.b*(xx2-fresult1.e))+fresult1.f;
    sig=[T1',T2];
    sig2=interp(sig,4); % plot(sig2)
    %     T1i=interp(T1,4); % RoR fit
    %     T2i=interp(T2,4); % tau fit
    %     bline2=min(sig); %%%%%%%%%%%% 191107
    %     AMP50=((sortida(6)-bline)/2)+bline2;
    %     AMP25=amp-((sortida(6)-bline)/4)+bline2;
    %     AMP75=amp-((sortida(6)-bline)*(3/4))+bline2;
    %     AMP90=amp-((sortida(6)-bline)*(90/100))+bline2;
    
    aux=find(sig2>AMP50); % >= en MVM
    if(~isempty(aux)),
        %plot(sig);
        dhmi=(aux(end)-aux(1))/4; % /(length(sig2)/length(sig)) --> Al interpolar tenim més punts, el DT és per la senyal original
        dhm=dhmi*DT;
        h50=AMP50;
        x51=x(1)+find(T(1:centre)<=AMP50,1,'last')*DT;
        x52=x51+round(dhm);
    else
        dhm=NaN;
    end
    % plot(sig2); hold on; scatter(aux(1),sig2(aux(1))); scatter(aux(end),sig2(aux(end)));
    
    sortida(5)=dhm;
    
    % Altres % de l'amp
    AMPP=[AMP25;AMP50;AMP75;AMP90];
    for kk=1:size(AMPP)
        AMPi=AMPP(kk);
        aux=find(sig2>=AMPi);
        if(~isempty(aux)),
            %plot(sig);
            dhmi=(aux(end)-aux(1))/4;
            dhmi=dhmi*DT;
            DHMi(kk)=dhmi;
            hi(kk)=AMPi;
            xpi(kk)=x(1)+find(T(1:centre)<AMPi,1,'last')*DT;
            xp2i(kk)=xpi(kk)+round(dhmi);
        else
            DHMi(kk)=NaN;
        end
    end
    sortida(7:10)=DHMi;
    
end




if(plt==1)
    
    y=Tfp;
    xlims=[min(xfp),max(xfp)];
    ylims=[-.1 normVal];
    
    marge=max([ceil(min(imsize)*.08),22]);
    I=zeros(imsize(1)-marge,imsize(2)-marge);
    if(strcmp(titol,'')~=1),I=zeros(imsize(1)-2*marge,imsize(2)-marge);end
    [a,b]=size(I);
    
    
    
    
    %calcula regio on pintar
    x0=[1:b];
    xR=linspace(xlims(1),xlims(2),b);
    if(xR(1)<x(1)),
        ori=find(xR>=x(1),1,'first');xR=xR(ori:end);x0=x0(ori:end);
    else
        ori=find(x>=xR(1),1,'first'); x=x(ori:end);y=y(ori:end);
    end
    if(xR(end)<x(end)),
        fin=find(x<=xR(end),1,'last');x=x(1:fin);y=y(1:fin);
    else
        fin=find(xR<=x(end),1,'last');xR=xR(1:fin);x0=x0(1:fin);
    end
    if(length(x0)>length(x)),
        x00=zeros(size(x));
        for ii=1:length(x)
            [aux,id]=min(abs(xR-x(ii)));
            x00(ii)=id;
        end
        x0=x0(x00);
        xR=xR(x00);
    else
        x00=zeros(size(x0));
        for ii=1:length(x0)
            [aux,id]=min(abs(x-xR(ii)));
            x00(ii)=id;
        end
        x=x(x00);
        y=y(x00);
    end
    
    
    yy=y;
    yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
    %yy(yy<=0)=1;yy(yy>=a)=a;
    % pinta grafic complet
    x00=1+round([0:length(xfp)-1]*(x0(end)-x0(1))/(length(x0)-1));
    for kk=2:min([length(x00),length(yy)])
        [allx1,ally1] = puntsEnmig([x00(kk-1) x00(kk)],[yy(kk-1) yy(kk)]);
        for kk2=1:length(allx1)
            lay=1+a-ally1(kk2);
            lax=allx1(kk2);
            if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b-gru+1)))
                I(lay,lax+gru-1)=5;
            end
        end
    end
    % pinta grafic bo
    for kk=2:length(x0)
        [allx1,ally1] = puntsEnmig([x0(kk-1) x0(kk)],[yy(kk-1) yy(kk)]);
        for kk2=1:length(allx1)
            lay=1+a-ally1(kk2);
            lax=allx1(kk2);
            if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b-gru+1)))
                I(lay,lax:lax+gru-1)=1;
            end
        end
    end
    %I=conv2(I,[1 1;1 1],'same');I(I>0)=1;
    
    cole=diff(ylims)*.01;% lu q sobresurten les linies horitzontals al final
    
    
    CP=0;col6=[];
    try
        yy=fitted2;
        rel=find((yy<=ylims(2))&(yy>=ylims(1)));
        yy=yy(rel);
        xx=x0(rel);
        yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
        %yy(yy<=0)=1;yy(yy>=a)=a;
        % pinta fit1 exp
        for kk=2:length(yy)
            [allx1,ally1] = puntsEnmig([xx(kk-1) xx(kk)],[yy(kk-1) yy(kk)]);
            for kk2=1:length(allx1)
                lay=1+a-ally1(kk2);
                lax=allx1(kk2);
                if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b)))
                    I(lay,lax)=2;
                end
            end
        end
        if(tau>1e3),tautx=num2str(tau,'%4.2E');else tautx=num2str(.01*round(tau*100));end
        CP=CP+1;parrafada{CP}=['tau = ' tautx 'ms (R²= ' num2str(.001*round(rsquare1*1000)) ')'];
        col6=[col6 2];
    end
    
    try
        yy=fitted1;
        rel=find((yy<ylims(2))&(yy>ylims(1)));
        yy=yy(rel);
        xx=x0(rel);
        yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
        %yy(yy<=0)=1;yy(yy>=a)=a;
        % pinta fit uprise
        for kk=2:length(yy)
            [allx1,ally1] = puntsEnmig([xx(kk-1) xx(kk)],[yy(kk-1) yy(kk)]);
            for kk2=1:length(allx1)
                lay=1+a-ally1(kk2);
                lax=allx1(kk2);
                if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b)))
                    I(lay,lax)=2;
                end
            end
        end
        CP=CP+1;parrafada{CP}=['ror = ' num2str(.001*round(ror*1000)) ' ms'];
        col6=[col6 2];
    end
    
    colex=abs(diff(xlims(1)+(abs((round(a*([0 cole]-ylims(1))/(ylims(2)-ylims(1)))))*(xlims(2)-xlims(1))/b)));
    rg=max(x)-min(x);ox=round(x(1)+.9*rg);
    
    try
        yy=[bline bline bline bline+amp bline+amp bline+amp];
        xx=[ox-colex ox+colex ox ox ox-colex ox+colex];
        xx=round(b*(xx-xlims(1))/(xlims(2)-xlims(1)));
        yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
        %yy(yy<=0)=1;yy(yy>=a)=a;
        % pinta linia vertical amplitud
        for kk=2:length(xx)
            [allx1,ally1] = puntsEnmig([xx(kk-1) xx(kk)],[yy(kk-1) yy(kk)]);
            for kk2=1:length(allx1)
                lay=1+a-ally1(kk2);
                lax=allx1(kk2);
                if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b)))
                    I(lay,lax)=4;
                end
            end
        end
        CP=CP+1;parrafada{CP}=['amp = ' num2str(.001*round(amp*1000))];
        col6=[col6 4];
    end
    try
        yy=[h50+cole h50-cole h50 h50 h50+cole h50-cole];
        xx=[x51 x51 x51 x52 x52 x52];
        xx=round(b*(xx-xlims(1))/(xlims(2)-xlims(1)));
        yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
        %yy(yy<=0)=1;yy(yy>=a)=a;
        % pinta linia fwhm
        for kk=2:length(xx)
            [allx1,ally1] = puntsEnmig([xx(kk-1) xx(kk)],[yy(kk-1) yy(kk)]);
            for kk2=1:length(allx1)
                lay=1+a-ally1(kk2);
                lax=allx1(kk2);
                if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b)))
                    I(lay,lax)=3;
                end
            end
        end
        CP=CP+1;parrafada{CP}=['FDHM = ' num2str(.01*round(dhm*100)) 'ms'];
        col6=[col6 3];
    end
    try
        pltd=1;
        if pltd==1,
            for kk=1:size(AMPP,1)
                yy=[hi(kk)+cole hi(kk)-cole hi(kk) hi(kk) hi(kk)+cole hi(kk)-cole];
                xx=[xpi(kk) xpi(kk) xpi(kk) xp2i(kk) xp2i(kk) xp2i(kk)];
                xx=round(b*(xx-xlims(1))/(xlims(2)-xlims(1)));
                yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
                %yy(yy<=0)=1;yy(yy>=a)=a;
                % pinta linia fwhm
                for kk=2:length(xx)
                    [allx1,ally1] = puntsEnmig([xx(kk-1) xx(kk)],[yy(kk-1) yy(kk)]);
                    for kk2=1:length(allx1)
                        lay=1+a-ally1(kk2);
                        lax=allx1(kk2);
                        if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b)))
                            I(lay,lax)=3;
                        end
                    end
                end
            end
        end
    end
    try
        yy=[bline+cole bline-cole bline bline bline+cole bline-cole];
        xx=[xb xb xb x(centre(1)) x(centre(1)) x(centre(1))];
        xx=round(b*(xx-xlims(1))/(xlims(2)-xlims(1)));
        yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
        %yy(yy<=0)=1;yy(yy>=a)=a;
        % pinta linia t2p
        for kk=2:length(xx)
            [allx1,ally1] = puntsEnmig([xx(kk-1) xx(kk)],[yy(kk-1) yy(kk)]);
            for kk2=1:length(allx1)
                lay=1+a-ally1(kk2);
                lax=allx1(kk2);
                if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b)))
                    I(lay,lax)=3;
                end
            end
        end
        CP=CP+1;parrafada{CP}=['t2p = ' num2str(.01*round(t2p*100)) 'ms'];
        col6=[col6 3];
    end
    
    
    % enxufa Y ticks
    [ytiks] = niceTicks(ylims);
    yimtiks=round(1+a-a*(ytiks-ylims(1))/(ylims(2)-ylims(1)));
    ytiks(yimtiks<1)=[];yimtiks(yimtiks<1)=[];
    ytiks(yimtiks>a)=[];yimtiks(yimtiks>a)=[];
    for jj=1:length(yimtiks)
        I(yimtiks(jj),1:round(marge*.5))=1;
    end
    
    
    % enxufa X ticks
    [xtiks] = niceTicks(xlims);
    ximtiks=round(b*(xtiks-xlims(1))/(xlims(2)-xlims(1)));
    xtiks(ximtiks<1)=[];ximtiks(ximtiks<1)=[];
    xtiks(ximtiks>b)=[];ximtiks(ximtiks>b)=[];
    for jj=1:length(ximtiks)
        I(a-round(marge*.5):a,ximtiks(jj))=1;
    end
    
    
    if(txtfunc~=1)% pas a rgb
        %temp=gaussiana2d(3);temp=temp/sum(sum(temp));
        I2=zeros(size(I));I2(:,:,2)=I2;I2(:,:,3)=zeros(size(I));
        %II=I;II=conv2(II,[1 1],'same');
        for ii=0:5
            for jj=1:3
                CC=I2(:,:,jj);
                CC(I==ii)=cols(ii+1,jj);
                I2(:,:,jj)=CC;%conv2(CC,[1 1],'same');
            end
        end
        I=I2;
    end
    % text dels parametres
    desf=0;
    for ii=1:CP
        if(txtfunc==1)
            Z=escriuFrase(parrafada{ii},FS);
            [zy,zx]=size(Z);
            ZB=I(3+1+desf:3+zy+desf,size(I,2)-zx+1-3:size(I,2)-3);
            Z(Z==0)=ZB(Z==0);
            Z(Z==1)=col6(ii);
            I(3+1+desf:3+zy+desf,size(I,2)-zx+1-3:size(I,2)-3)=Z;
            desf=desf+zy;
        else
            I=textIm(size(I,2)+1-3,3+1+round(FS*.8*(ii-1)),parrafada{ii},double(I),'fontsize',FS,'textcolor',cols(col6(ii)+1,:),'horizontalalignment','right','verticalalignment','top','blending','on');
        end
    end
    
    % marc
    if(txtfunc==1)
        I(1,:)=1;I(end,:)=1;I(:,1)=1;I(:,end)=1;
    else
        for ii=1:3
            I(1,:,ii)=cols(2,ii);I(end,:,ii)=cols(2,ii);I(:,1,ii)=cols(2,ii);I(:,end,ii)=cols(2,ii);
        end
    end
    % afegeix marge
    if(txtfunc==1)
        I=[zeros(a,marge) I];I=[I;zeros(marge,b+marge)];
    else
        I=cat(2,ones(a,marge,3),I);I=cat(1,I,ones(marge,b+marge,3));
    end
    % enxufa text per Y tiks
    for jj=1:length(yimtiks)
        if(txtfunc==1)
            Z=escriuFrase(num2str(ytiks(jj)),FS);
            if(size(Z,2)<marge)
                oriy=round(yimtiks(jj)-size(Z,1)/2);
                finy=oriy+size(Z,1)-1;
                orix=marge-size(Z,2);
                finx=orix+size(Z,2)-1;
                if(oriy>0)
                    I(oriy:finy,orix:finx)=Z;
                end
            end
        else
            I=textIm(marge,yimtiks(jj),num2str(ytiks(jj)),I,'fontsize',FS,'textcolor',cols(2,:),'horizontalalignment','right','verticalalignment','middle','blending','on');
        end
    end
    
    % enxufa text per X tiks
    for jj=1:length(ximtiks)
        if(txtfunc==1)
            Z=escriuFrase(num2str([num2str(xtiks(jj))]),FS);
            orix=round(marge+ximtiks(jj)-size(Z,2)/2);
            finx=orix+size(Z,2)-1;
            oriy=a+1;
            finy=oriy+size(Z,1)-1;
            if(finx<b+marge)
                I(oriy:finy,orix:finx)=Z;
            end
        else
            I=textIm(marge+ximtiks(jj),a+1,num2str(xtiks(jj)),I,'fontsize',FS,'textcolor',cols(2,:),'horizontalalignment','center','verticalalignment','top','blending','on');
        end
    end
    % enxufa titol si cal
    if(strcmp(titol,'')~=1)
        if(txtfunc==1)
            I=[zeros(marge,b+marge);I];
            Z=escriuFrase(titol,FS);
            oriy=round((marge/2)-(size(Z,1)/2));
            orix=round((size(I,2)/2)-(size(Z,2)/2));
            I(oriy:oriy+size(Z,1)-1,orix:orix+size(Z,2)-1)=Z;
        else
            I=cat(1,ones(marge,b+marge,3),I);
            I=textIm(round(size(I,2)/2),round(marge/2),titol,I,'fontsize',FS+2,'textcolor',cols(2,:),'horizontalalignment','center','verticalalignment','middle','blending','on');
        end
    end
    % xlabel
    if(txtfunc==1)
        Z=escriuFrase('time (ms)',FS);
        oriy=round(size(I,1)-(size(Z,1)));
        orix=round((size(I,2)/2)-(size(Z,2)/2));
        I(oriy:oriy+size(Z,1)-1,orix:orix+size(Z,2)-1)=Z;
    else
        I=textIm(round(size(I,2)/2),size(I,1),'time (ms)',I,'fontsize',FS,'textcolor',cols(2,:),'horizontalalignment','center','verticalalignment','bottom','blending','on');
    end
    % % ylabel
    % if(strcmp(ylabl,'')~=1)
    %  Z=escriuFrase(ylabl,FS);
    %  Z=Z';Z=Z(end:-1:1,:);
    % orix=2;
    % oriy=round((size(I,1)/2)-(size(Z,2)/1));
    % I(oriy:oriy+size(Z,1)-1,orix:orix+size(Z,2)-1)=Z;
    % end
    
    % espai al final
    if(txtfunc==1)
        I=[I zeros(size(I,1),1)];
    else
        I=cat(2,I,ones(size(I,1),1,3));
    end
    
    
    
end






end



% AMP
%     [amp ind]=max(T(centre(1)-3:centre(1)+8));% (F-Fo)/Fo
%     centre(1)=centre(1)-3+ind-1; %%%%%%%%% CORRECCIO MAX!
%     bline=quantile(T(1:centre(1)),.02);
%     sortida(6)=amp;
%     amp=amp-bline;
%     sortida(1)=amp;

%%%%%% RoR
%     % tram de pujada
%     aux1=centre(1);
%     if(aux1>=10),aux1=find(T>T(centre(1))*.9,1,'first');end
%     uprise=T(1:aux1);
%     di=zeros(1,aux1-1);
%     for ii=1:aux1-1
%         t1=uprise(1:ii);
%         t2=uprise(ii:centre(1));
%         [fresult2,S]=polyfit([ii:centre(1)],t2,1);
%         t3=[ii:centre(1)]*fresult2(1)+fresult2(2);
%         %     figure;plot(uprise);hold on;
%         %     plot([ii:centre(1)],t3,'r');
%         %     plot([1:ii],mean(t1)*ones(1,ii),'r');
%         ft=[mean(t1)*ones(1,ii) [ii:aux1]*fresult2(1)+fresult2(2)];
%         ft(ii)=(ft(ii)+ft(ii+1))/2;
%         ft(ii+1)=[];
%                di(ii)=abs(mean(t1)-t3(1));% distancia entre punts dencaix dels fits
%         % di(ii)=sum(abs(ft'-uprise)); % distancia entre senyal i fits
%     end
%     [kk,aux]=min(di);
%
%     T1=T(aux:aux1);
%     x1=x(aux:aux1);
%     %figure;plot(x,T);hold on;plot(x1,T1,'r');
%
%     % fit lineal a la pujada
%     if(size(x1,1)<size(x1,2)),x1=x1';end
%     if(size(T1,1)<size(T1,2)),T1=T1';end
%     [fresult2,S]=polyfit(x1,T1,1);
%     T1=fresult2(1)*x1+fresult2(2);
%     fitted1=fresult2(1)*x+fresult2(2);
%     %figure;plot(x,fitted1);hold on;plot(x,fitted2);
% RoR
%     % tram de pujada
%     % %%%% abans
%     % uprise=T(1:centre(1));
%     % du=diff(uprise);
%     % [aux dd]=max(du);
%     % du=smooth(du,round(length(du)/8));
%     % lix=zeros(length(du),1);
%     % lix(du>max(du)*.3)=1;
%     % liz=1-lix;liz=bwlabel(liz);liz2=zeros(size(liz));
%     % for ii=1:max(liz),
%     % if(numel(find(liz==ii))<=3),liz2(liz==ii)=1;end;
%     % if(find(liz==ii,1,'first')>length(liz)*.7),liz2(liz==ii)=1;end
%     % end
%     % lix(liz2==1)=1;lix=bwlabel(lix);
%     % lix(lix~=lix(find(lix>0,1,'last')))=0;lix(lix>0)=1;
%     % lix=find(lix>0);
%     % aux=lix(1);%aux=lix(ceil(length(lix)*.1))+1;
%     % aux1=lix(end)+1;%aux1=lix(ceil(length(lix)*.9))+1;
% %%%%% ara
%     aux1=centre(1);
%     uprise=T(1:aux1);
%     di=zeros(1,aux1-1);
%     for ii=1:aux1-1
%         t1=uprise(1:ii);
%         t2=uprise(ii:centre(1));
%         [fresult2,S]=polyfit([ii:centre(1)],t2,1); % [fresult2,S]=polyfit([ii:centre(1)]',t2,1);
%         t3=[ii:centre(1)]*fresult2(1)+fresult2(2);
%         %     figure;plot(uprise);hold on;
%         %     plot([ii:centre(1)],t3,'r');
%         %     plot([1:ii],mean(t1)*ones(1,ii),'r');
%         ft=[mean(t1)*ones(1,ii) [ii:centre(1)]*fresult2(1)+fresult2(2)];
%         ft(ii)=(ft(ii)+ft(ii+1))/2;
%         ft(ii+1)=[];
%
%         di(ii)=sum(abs(ft-uprise)); % di(ii)=sum(abs(ft'-uprise));
%     end
%     [kk,aux]=min(di);
%%%%%% FDHM
%     aux=centre(1)-1+find(T(centre(1):end)<bline+amp/2,1,'first');
%     aux2=find(T(1:centre(1))<bline+amp/2,1,'last');
%     AMP50=(bline+amp)/2;
%     if((~isempty(aux))&&(~isempty(aux2))),dhm=x(aux)-x(aux2);x51=x(aux2);x52=x(aux);h50=T(aux); else dhm=NaN;end
% %     %  CARMETA 181011
% %      if((~isempty(aux))&&(~isempty(aux2))),
% %          if T(aux) ~= AMP50, aux3=aux+1;  m=(T(aux3)-T(aux))/(aux3-aux); auxf=((AMP50-T(aux))/m)+aux3; else auxf=aux;end
% %          if T(aux2) ~= AMP50, aux4=aux2-1;  m=(T(aux4)-T(aux2))/(aux4-aux2); auxi=((AMP50-T(aux2))/m)+aux4; else auxi=aux2;end
% %          dhm=(auxf-auxi)*DT; x51=x(1)+auxi*DT;x52=x(1)+auxf*DT;h50=AMP50; 
% %          aux=auxf;aux2=auxi;
% %      else
% %          dhm=NaN;
% %      end
% %     
% %     %  CARMETA

    % aux=centre(1)-1+find(T(centre(1):end)<bline+amp/2,1,'first');
%     aux2=find(T(1:centre(1))<bline+amp/2,1,'last');
%     AMP50=(bline+amp)/2;
%%%%%% tau
% fit exponencial al decay
%     model1 = fittype('a*exp(-b*(x-e))+f');
%     options1 = fitoptions(model1);
%     options1.Display = 'off';
%     options1.Robust = 'bisquare';
%     %%%%%%%%%%%%%%%%%%%%%%      a         b        e        f
%     options1.startpoint = [max(abs(TT)) 0  min(xx) min(TT)];
%     options1.Lower = [0 0 -10*max(xx) min(TT)];
%     options1.Upper = [100*max(abs(TT)) 1 100*max(xx) max((TT))];
%     [fresult1,gof2,out2]  = fit(xx',TT,model1,options1);
%     rsquare1=gof2.rsquare;
%     r2=rsquare1;
%     %     figure;plot(x,T,'-b'); hold on;plot(fresult1,'r');
%     fitted2=fresult1.a*exp(-fresult1.b*(x-fresult1.e))+fresult1.f;
%     tau=1/fresult1.b;
%     sortida(2)=tau;
%     %sortida=[amp,tau,ror,t2p,dhm];




