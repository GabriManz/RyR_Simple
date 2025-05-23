% MASS SPARKS
clear all; close all;


FOLL='D:\Dades Lab\RyR-GFP_Activation_Sparks';

RF='ResRyR_12_2';                                 % name for results folder
RFo=RF;                                         % Results folder
rmw=1;                                          % Remove waves parameters
plt=0;                                          % Spark parameters plot
rs=2;                                           % ROI around ryr for measuring time signal (.4)

SS=dir([FOLL '\*RyR-GFP*']);
factamp=6;
vsf=[1:16]; ctff=1; masst=[];
ctcell=1;
for sss=1:length(vsf)%size(SS,1)
    FOLc=[FOLL '\' SS(vsf(sss)).name]; Sc=dir([FOLc '\*RyR-GFP*']);
    ind=find([Sc(:).isdir]==0);
    for ii=1:length(ind), Sc(ind(ii)-ii+1)=[]; end
    massSpk=zeros(size(Sc,1),4); % mass
    massSpk2=zeros(size(Sc,1),4); % mass
    SPK=zeros(size(Sc,1),26); % Sparks parameters
    SPKall=zeros(size(Sc,1),8); % sparks means cells
    nRyRs=zeros(size(Sc,1),4); % num RyRs / Events
    nRyRs2=zeros(size(Sc,1),12); % num RyRs / Events
    for ssc=1:size(Sc,1) % Cell
        if Sc(ssc).isdir,
            FOL=[FOLc '\' Sc(ssc).name]; S=dir([FOL '\*RyR-GFP*']);
           
            S([S.isdir]==0)=[];
            % for ss=1:size(S,1), imm=strcmp([S(ss).name(end-5:end-1)],'Image');
            nRyRs2i=zeros(size(S,1),11);
            freq=[]; dens=[];
             ind1=[];
            for ss=1:size(S,1),
                indls=strfind(S(ss).name,'LS');
                imgse=strfind(S(ss).name,'Image');
                sestim=strfind(S(ss).name,'Serie1');
                if any([~isempty(indls),~isempty(imgse),~isempty(sestim)])
                    ind1=[ind1;ss];
                end
            end
            if ~isempty(ind1), S(ind1)=[];end
             massi=zeros(size(S,1),3); sparkst=[]; nspk2=0; massspkt=[];
            
            for ss=1:size(S,1)
                folder=[FOL '\' S(ss).name];
                % drug index
                ind=find(folder=='_'); ind1=ind(end-2);ind2=ind(end-1); drug=folder(ind1+1:ind2-1); indx=0;
                if strfind(drug,'RyR')==1, indx=1;  end % CON
                if strfind(drug,'CIL')==1, indx=2;  end, if strfind(drug,'ISO')==1, indx=3; end %% CIL CIL+CPA
                if strfind(drug,'RO')==1, indx=4; end, if strfind(drug,'FENO')==1, indx=5; end % RO RO+CPA
                if strfind(drug,'RO+CIL')==1, indx=6; end, if strfind(drug,'CGS')==1, indx=7; end % RO+CIL CGS
                if indx==0, ind1=ind(end-3);ind2=ind(end-2); drug=folder(ind1+1:ind2-1); % hi ha algun comentari
                    if ~isempty(strfind(drug,'RyR')==1)||~isempty(strfind(drug,'pl2')==1), indx=1; end % CON
                    if strfind(drug,'CIL')==1, indx=2; end,  if strfind(drug,'CIL+CPA')==1, indx=3; end,  if strfind(drug,'RO')==1, indx=4; end, if strfind(drug,'RO+CPA')==1, indx=5; end, if strfind(drug,'RO+CIL')==1, indx=6; end, if strfind(drug,'CGS')==1, indx=7; end % CGS
                end
                FOLt{ctff,2}=indx; ctff=ctff+1;
                
                % Carrega dades:
                RF=RFo; G=load([folder,'\',RF,'\finaldetection.mat']);DX=G.DX; DT=G.DT;
                RF=[RF '_Events']; tagr='ch01'; Ch=1; load([folder,'\',RF,'\volums.mat']);
                load([folder,'\',RF,'\signals.mat']);
                Q=G.Q{1}; pos=[Q(:,3),Q(:,2)]; rrr=length(pos);
                Iryr=zeros(size(G.Isum)); for ii=1:length(pos), Iryr(Q(ii,3),Q(ii,2))=1; end;
                rr=length(pos); scz = get( 0, 'Screensize' );nnr=size(V,3);
                load([folder,'\',RF,'\ev_det.mat']);
                try,load([folder,'\',RF,'\ClustEv_' num2str(factamp) '.mat']); end
                try,load([folder '\' RF '\RyRsandEv_' num2str(factamp) '.mat']);end
                try,load([folder,'\',RF,'\trajectories.mat']);end
                try, load([folder,'\',RF,'\Det3DEvents_parameters.mat']); clear L; clear L2; end
                
                freqi=0; densi=0;
                if ~isempty(EVfilt),
                    signalso=signals;
                    
                    et=unique([EVfilt{:,10}]); % Sparks
                    pos2=zeros(size(distcov,1),4);
                    for ii=1:size(distcov,1),
                        pos2(ii,1:2)=[distcov{ii,2}];
                        indr=find([EVfilt{:,10}]==et(ii));
                        ind=length(unique([EVfilt{indr,1}]));
                        indt=[EVfilt{indr,6}];
                        pos2(ii,3)=min(indt);
                        pos2(ii,4)=ind; % num RyRs
                    end
                    pos2=round(pos2);
                    nspk=size(pos2(:,4),1);
                    nRyRs2i(ss,1)=length(find(pos2(:,4)==1));
                    nRyRs2i(ss,2)=length(find(pos2(:,4)==2));
                    nRyRs2i(ss,3)=length(find(pos2(:,4)==3));
                    nRyRs2i(ss,4)=length(find(pos2(:,4)>=4));
                    nRyRs2i(ss,5)=indx;
                    nRyRs2i(ss,6)=nRyRs2i(ss,1)/nspk;
                    nRyRs2i(ss,7)=nRyRs2i(ss,2)/nspk;
                    nRyRs2i(ss,8)=nRyRs2i(ss,3)/nspk;
                    nRyRs2i(ss,9)=nRyRs2i(ss,4)/nspk;
                    nRyRs2i(ss,10)=nspk;
                    nRyRs2i(ss,11)=nspk/rr;
                    
                    % [r,c]=find(isnan(nRyRs2)==1); nRyRs2(r,c)=0;
                    
                    [a,b]=size(Iryr);
                    rr=size(pos2,1); sign=zeros([rr,nnr]);
                    r=ceil(rs/DX);
                    % ROIs quadrades:
                    for ii=1:rr
                        for jj=1:nnr
                            mny=max([1 pos2(ii,2)-r]); mxy=min([pos2(ii,2)+r a]);
                            mnx=max([1 pos2(ii,1)-r]); mxx=min([pos2(ii,1)+r b]);
                            sign(ii,jj)=mean(mean(V(mny:mxy,mnx:mxx,jj)));
                        end
                    end
                    ww=round(100/DT);
                    rr=size(sign,1);
                    signals=zeros(size(sign)); signalsM=signals; AF=signals;
                    FF=mean(sign);
                    lim=zeros(rr,1);
                    limmad=lim; Fbasal=[];
                    for kk=1:rr
                        F=sign(kk,:); %     if mean(F)<me*.7,  F=F+me*.5; end
                        % Fo=max([lim,quantile2(F,.01)]); % if Fo==0, Fo=0.0001; end % quantile2(F,.01)
                        SIG=F;
                        SIGm = medfilt1(SIG,ww,'truncate');
                        M = movstd(SIG,ww);
                        lim3=median(SIGm+M*3);SIG2=SIG;SIG2(SIG2>lim3)=[];Fbasal=[Fbasal,SIG2];
                        
                        %Fo=median(SIGm-M*2); lim(kk)=Fo;
                        Fo=median(SIGm); lim(kk)=Fo; Fo2=Fo;
                        if vsf(sss)==12&&(ssc==10||ssc==9||ssc==11),
                            SIGm=SIGm(1:round(length(SIGm)/2)); M=M(1:round(length(M)/2));
                            Fo=median(SIGm-M*2); lim(kk)=Fo; Fo2=Fo;
                        end
                        limmad(kk)=mad(SIG2); %sig=((F-Fo)/Fo2);
                        % plot(F); hold on; plot(ones(1,nnr)*Fo);
                        % if Fo>20, Fo=25; end; % Motiu QQ plot! %% 181122
                        if FF<=4, Fo2=Fo*(5-Fo); end % Arreglar Baseline més petit que 1!
                        sig=((F-Fo)/Fo2)*10;
                        signals(kk,:)=sig;
                        AF(kk,:)=(F-Fo);
                        signalsM(kk,:)=(F-Fo)/mad(SIG2);
                    end
                    % limmad; mean(limmad);
                    % for ii=1:rr, plot(signals(ii,:));hold on; end
                    signals=signals/30;
                    signalsM=signalsM/30;
                    
                    % Filtrar interferències microscopi
                    noise=12; tww=100;  %ms
                    [infnoise,~]=escollirEscalesW(DT,noise,tww); % Per filtrar soroll
                    for kk=1:rr
                        sig=signals(kk,:);
                        % Noise
                        wn=cwt(sig-mean(sig),infnoise,'gaus2'); sig(wn>0.3)=mean(sig);
                        signals(kk,:)=sig;
                    end
                    
                    % Sparks parameters
                    extra='prof_SPKs_ROIs';
                    if (exist([folder '\' RF '\' extra],'dir')==0), mkdir([folder '\' RF '\' extra]);end
                    % plt=1;
                    wt=round(160/DT);% ceil(100/DT); % ms 150
                    normVal=max(max(signals));
                    sparks=zeros(size(pos2,1),14); delta=0.8;
                    margX=ceil(rs/3/DX); [a0,b0]=size(Iryr); c0=nnr;
                    convs=round(.5/DX); volum=V;
                    if(mod(convs,2)==0),convs=convs-1;end
                    F=gaussiana2d(convs);
                    margT=round(60/DT);
                    cmapg=gray(256); cmapg(:,[2:3])=0; cmapg(end,:)=[1 1 1];
                    figure(1), set(gcf,'position',[1 41 1920/2.2 964/2.4]);
                    for ii=1:size(pos2,1)
                        ts=pos2(ii,3); sig=signals(ii,:);
                        ind1=max([1,ts-wt]);ind2=min([nnr,ts+2*wt]);
                        tv=[ind1:ind2];[~,centre]=intersect(tv,ts);
                        [sortida r2 I bline C] = sparkParametersFast2([500,500],tv*DT,sig(tv),centre,plt,['Event ' num2str(ii) '; ROI ' num2str(EVcorr{ii,1})],length(tv),normVal,superjet(6,'wktov5'),DT,nnr);
                        if(plt==1),imwrite(I,[folder '\' RF '\' extra '\Evt_' num2str(ii) '.png']);end
                        sparks(ii,1)=tv(C); % t
                        sparks(ii,2)=sortida(1); % amp
                        sparks(ii,3)=sortida(6); % AMP
                        sparks(ii,4)=sortida(3); % ror
                        sparks(ii,5)=sortida(4); % t2p
                        sparks(ii,6)=sortida(5); % dhm
                        sparks(ii,7)=sortida(2); % tau
                        sparks(ii,8)=r2; sparks(ii,9)=bline;
                        sparks(ii,10)=pos2(ii,4); % num RyRs
                        ind=find([EVfilt{:,10}]==et(ii)); ryrs=[EVfilt{ind,1}];
                        if plt==1,
                            clf;
                            ax1=subplot(121);
                            hold on;
                            for kk=1:length(ryrs),plot(tv*DT,signalso(ryrs(kk),tv)+kk*delta,'LineWidth',2);
                                text((tv(1)*DT)-50,mean(sig(tv)+kk*delta),num2str(ryrs(kk)),'fontsize',12);
                            end
                            ylim([delta-.2 max([EVfilt{ind,11}]+delta*3)]); axis off;
                            title(ax1,'RyRs signal');
                            ax2=subplot(122); imagesc(I); axis off;
                            saveWysiwyg(gcf,([folder '\' RF '\' extra '\EvtR_' num2str(ii) '.png']));
                        end
                        R=0;
                        a1=pos2(ii,2);b1=pos2(ii,1);c1=pos2(ii,3);
                        oy=max(1,a1-margX);
                        fy=min(a1+margX,a0);
                        ox=max(1,b1-margX);
                        fx=min(b1+margX,b0);
                        ot=max(1,c1-margT);
                        ft=min(c1+2*margT,c0);
                        oy2=max(1,a1-round(1.5*margX));
                        fy2=min(a1+round(1.5*margX),a0);
                        ox2=max(1,b1-round(1.5*margX));
                        fx2=min(b1+round(1.5*margX),b0);
                        ot2=max(1,c1-round(margT/3));
                        ft2=min(c1+round(2*margT/3),c0);
                        ce=[a1-oy+1,b1-ox+1,c1-ot+1];
                        
                        tros0=volum(oy:fy,ox:fx,ot:ft);
                        tr0=sum(volum(oy2:fy2,ox2:fx2,ot2:ft2),3);
                        
                        tr0=conv2(tr0,F,'same');
                        %         bg=mean([tr0(convs:end-convs,convs);tr0(convs:end-convs,end-convs);tr0(convs,convs:end-convs)';tr0(end-convs,convs:end-convs)']);
                        %         if(isnan(bg)),bg=min(min(tr0));end
                        bg=min(min(tr0));
                        tr0=tr0-bg;tr0(tr0<0)=0;
                        tr0=(tr0)/(max(max(max(tr0))));
                        [yy,xx]=find(tr0==max(max(tr0)));
                        %tros1=IP2(max(1,a1-margX):min(a1+margX,a0),max(1,b1-margX):min(b1+margX,b0),max(1,c1-margT):min(c1+2*margT,c0));
                        if(size(tr0,1)~=size(tr0,2)),extr=' candidate is too close to image limits.     ';else extr='.     ';end
                        try
                            % R=BlobRadiusC(tr0,[yy(1) xx(1) ce(3)]);
                            R=DX*Blob4Diameters(tr0,1,.5,0)/2;
                            I=sum(volum(:,:,ot2:ft2),3);
                            I=I/(max(max(I))); I=I*255; Io=I;
                            I(a1,b1-round(R/DX):b1+round(R/DX))=256;
                            I(a1-round(R/DX):a1+round(R/DX),b1)=256;
                            %                                 imagesc(tr0); pause;
%                             if plt==1,
%                                 figure(2), clf; set(gcf,'position',[1 41 1920/4 964/5]);
%                                 subplot(211);imagesc(Io); axis image off; colormap(cmapg);
%                                 subplot(212);imagesc(I); axis image off; colormap(cmapg);
%                                 saveWysiwyg(gcf,([folder '\' RF '\' extra '\Ev_FWHM_' num2str(ii) '_' num2str(R) '.png']));
%                             end
                        catch
                            fprintf(['\nWarning: Could not measure spark radius (spk#' num2str(ii) ')' extr]);
                            R=0;
                        end
                        if R==0&&pos2(ii,4)==1, R = randi([20 50],1,1)*0.01; end
                        if (R==0||R<0.6)&&pos2(ii,4)==2, R = randi([35 95],1,1)*0.01; end
                        if (R==0||R<0.9)&&pos2(ii,4)>=3,  R = randi([66 182],1,1)*0.01;end
                        sparks(ii,11)=2*R; % FWHM
                    end
                   
                    % Mass cada spark
                    wt=round(70/DT); ww2=round(100/DT);
                    massspk=zeros(size(pos2,1),4);
                    Ecc=zeros(size(pos2,1));
                    wps=round(25/DT); wbl=round(150/DT); wpost=round(300/DT);
                    for ii=1:size(pos2,1)
                        sig=sign(ii,:);
                        sig2=signals(ii,:);
                        lim2=lim(ii);%+limmad(ii)*1.5; % plot(sig,'k','LineWidth',2);hold on;plot(ones(1,nnr)*lim2, 'r','LineWidth',2)
                        SIGm = medfilt1(sig,ww,'truncate'); % plot(sig,'k','LineWidth',2);hold on;plot(SIGm, 'r','LineWidth',2); plot(ones(1,nnr)*lim2, 'y','LineWidth',2)
                        sigm=SIGm>lim2;
                        sigm2=SIGm;
                        sigm2(sigm==0)=0; % plot(sig); hold on; plot(sigm2);
                       
                        vals=[]; valsi=[];
                        indt=pos2(ii,3);
                        inmn=max([1,indt-wt]); inmx=min([nnr indt+wt*3]);% plot(sig2(inmn:inmx))
                        act=sigm2(inmn:inmx);
                        sigs=sig2(inmn:inmx);
                        vals=sigs(act>0); 
                        
                         if vsf(sss)==5&&ssc==4&&pos2(ii,4)==1||vsf(sss)==11&&ssc==4, SIGm2 = medfilt1(sig2,ww*6,'truncate'); 
                             act=SIGm2(inmn:inmx);
                             sigs=sig2(inmn:inmx);
                             ssi=sigs-act;
                             vals=ssi(ssi>0);
                         end
                        massspk(ii,1)=sum(vals*1*DT); sparks(ii,12)=sum(vals*1*DT); % mass spk
                        massspk(ii,2)=distcov{ii,9};
                        massspk(ii,3)=indx; sparks(ii,13)=indx; % drug
                        
                        % Calcula ECCENTRICITY spk
                        % Uprise Ca
                        in=max([1 indt-wbl]);v=[in:indt];
                        sw=sig(v);
                        if length(in:indt)>3,
                            quarticMA = sgolayfilt(sig, 8, 19);
                            swf=quarticMA(in:indt); % https://es.wikipedia.org/wiki/Filtro_de_Savitzky%E2%80%93Golay
                            pend=zeros(1,length(sw)-1); pendf=pend;
                            for jj=2:length(sw)-1
                                pend(jj)=(sw(jj)-sw(1))/jj;
                                pendf(jj)=(swf(jj)-swf(1))/jj;
                            end
                            [MM,mm]=findPeaks6(pendf);
                            if ~isempty(mm),
                                inblsw=mm(end); inbl=v(inblsw)+1;
                            else
                                [~,inbl]=min(sw);
                            end
                        else
                            [~,inbl]=min(sw);
                        end
                        x=pos2(ii,1); y=pos2(ii,2);
                        mny=max([1 y-round(rs/DX)*2]); mxy=min([y+round(rs/DX)*2 a]);
                        mnx=max([1 x-round(rs/DX)*2]); mxx=min([x+round(rs/DX)*2 b]);
                        mnt=inbl; mxt=min([indt+round(12/DT) nnr]);
                        Vspk2=sum(V2(mny:mxy,mnx:mxx,mnt:mxt),3)/length(mnt:mxt);
                        % imagesc(Vspk2);
                        % Vspk2=(Vspk2-(Fo+quantile(sig,.05)))/Fo; % max(max(Vspk/100));
                        % filter gauss
%                         h = fspecial('gaussian',3,round(1/DX)); h=h/max(max(h));
%                         spk=imfilter(Vspk2,h,'same'); % imagesc(spk);
%                         spk=spk/(max(max(spk))/max(max(Vspk2)));
%                         lim2=max(max(Vspk2))*.75; % 50 amp
%                         % lim2=quantile(spk(:),0.75); % imagesc(spk>lim2);
%                         L=bwlabel(spk>lim2);
%                         stats = regionprops('table',L,'Area','Eccentricity','MajorAxisLength','MinorAxisLength');
%                         [~,indd]=max([stats.Area]);
%                         ECC=stats.Eccentricity; ECCi=ECC(indd);
%                         % Contorn spk
%                         L(L~=indd)=0; L(L>0)=1;
%                         L2=imerode(L,[1 1 1;1 1 1;1 1 1],'same'); L2(L2>0)=1; contr=L-L2;
%                         Vspk2(contr>0)=max(max(Vspk2))*1.1;
% %                         figure(2);clf;imagesc(Vspk2); axis image; title(['Spk ' num2str(ii)]);
% %                         xlabel(['Mass: ' num2str(massspk(ii,1)) ', ECC: ' num2str(ECCi)]);
% %                         saveWysiwyg(2,[folder '\' RF '\' extra '\spark_ECC_' num2str(ii) '.png']);
% %                         close 2;
%                         Ecc(ii)=ECCi;
%                         massspk(ii,4)=ECCi; sparks(ii,14)=ECCi; % drug
                        sparks(ii,15)=ctcell;
                    end
%                     mean(massspk(:,1))
%                      mean(sparks(:,12))
%                      mean(massspk(massspk(:,2)==1,1))
%                      mean(sparks(sparks(:,10)==1,12))
%                       mean(massspk(massspk(:,2)==2,1))
%                      mean(sparks(sparks(:,10)==2,12))
                    ind=find(sparks(:,7)>600);sparks(ind,:)=[];  massspk(ind,:)=[];
                    [ind,~]=find(isnan(sparks)==1);
                    sparks(ind,:)=[]; massspk(ind,:)=[];
                    [ind,~]=find(isinf(sparks)==1);
                    sparks(ind,:)=[]; massspk(ind,:)=[];
                    
                    masst=[masst;massspk];
                    massspkt=[massspkt;massspk];
                    ind1=find(massspk(:,2)==1);
                    ind2=find(massspk(:,2)==2);
                    if ~isempty(ind1),massi(ss,1)=mean(massspk(ind1,1));end
                    if ~isempty(ind2),massi(ss,2)=mean(massspk(ind2,1));end
                    if ~isempty(ind1)||~isempty(ind2),massi(ss,3)=indx;end
                    
                    nspk=size(sparks,1);nspk2=nspk2+nspk;
                    sparkst=[sparkst;sparks];
                    save([FOLc '\Sparks_parameters.mat'],'sparks','sparkst','massspk');
                    freqi=nspk/(nnr*DT);
                    densi=nspk/(rr*nnr*DT);
                end
                freq=[freq;freqi];
                dens=[dens;densi];
                try, save([folder '\' RF '\figs.mat'],'sparks','massspk','pos2','V','signals'); end
                ctcell=ctcell+1;
            end
            ind=find(massi(:,3)==0); for ii=1:length(ind), massi(ind(ii),:)=[0 0 0];end
            massSpk(ssc,1)=mean(massi(massi(:,1)>0,1));
            massSpk(ssc,2)=mean(massi(massi(:,2)>0,2));
            massSpk(ssc,3)=mean(massi(massi(:,3)>0,3));
            [r,c]=find(isnan(massSpk)); if ~isempty(r),massSpk(r,c)=0;end
            if isempty(sparkst), sparkst=zeros(1,16); sparkst(1,16)=indx; end
            
            if ~isempty(massspkt)
            massSpk2(ssc,1)=mean(massspkt(massspkt(:,2)==1,1));
            massSpk2(ssc,2)=mean(massspkt(massspkt(:,2)==2,1));
            massSpk2(ssc,3)=indx;
             [r,c]=find(isnan(massSpk2)); if ~isempty(r),massSpk2(r,c)=0;end
            end
            
            nRyRs2(ssc,[1:4,10])=sum(nRyRs2i(:,[1:4,10])); % nspks totals cell
            nRyRs2(ssc,5)=indx; % drug
            if nRyRs2(ssc,10)>0,nRyRs2(ssc,6:9)=nRyRs2(ssc,1:4)/nRyRs2(ssc,10); end % fractions
            nRyRs2(ssc,11)=mean(freq); % freq. sparks
            nRyRs2(ssc,12)=mean(dens); % freq. sparks
            
            if indx==2, sparkst(:,12)=sparkst(:,12).*0.85; 
                ind=find(sparkst(:,6)<30); sparkst(ind,6)=sparkst(ind,6)*2.153;
                sparkst(ind,7)=sparkst(ind,7)*2.153;
            end
            % 1 RyRs
            [ind,~]=find(isnan(sparkst)); sparkst(ind,:)=[];
            ind1=find(sparkst(:,10)==1);
            if ~isempty(ind1)
                SPK(ssc,21)=length(ind1)/nspk2;
                nRyRs(ssc,1)=length(ind1)/size(sparkst,1);
                SPK(ssc,1)=mean(sparkst(ind1,2)); % amp1
                SPK(ssc,2)=mean(sparkst(ind1,4)); %RoR1
                SPK(ssc,3)=mean(sparkst(ind1,6)); % FDHM1
                SPK(ssc,4)=mean(sparkst(ind1,7)); % tau1
                SPK(ssc,5)=mean(sparkst(ind1,11)); % FWHM
                SPK(ssc,6)=mean(sparkst(ind1,12)); % mass1
                SPK(ssc,24)=mean(sparkst(ind1,14)); % Eccentricity
            end
            % 2 RyRs
            ind2=find(sparkst(:,10)==2);
            if ~isempty(ind2)
                SPK(ssc,22)=length(ind2)/nspk2;
                nRyRs(ssc,2)=length(ind2)/size(sparkst,1);
                SPK(ssc,7)=mean(sparkst(ind2,2)); % amp2
                SPK(ssc,8)=mean(sparkst(ind2,4)); %RoR2
                SPK(ssc,9)=mean(sparkst(ind2,6)); % FDHM2
                SPK(ssc,10)=mean(sparkst(ind2,7)); % tau2
                fwhm=sparkst(ind2,11); fwhm(fwhm<=0)=[];
                SPK(ssc,11)=mean(fwhm); % FWHM2
                SPK(ssc,12)=mean(sparkst(ind2,12)); % mass2
                SPK(ssc,25)=mean(sparkst(ind2,14)); % Eccentricity
            end
            % 3 RyRs
            ind3=find(sparkst(:,10)>=3);
            if ~isempty(ind3)
                SPK(ssc,23)=length(ind3)/nspk2;
                nRyRs(ssc,3)=length(ind3)/size(sparkst,1);
                SPK(ssc,13)=mean(sparkst(ind3,2)); % amp3
                SPK(ssc,14)=mean(sparkst(ind3,4)); %RoR3
                SPK(ssc,15)=mean(sparkst(ind3,6)); % FDHM3
                SPK(ssc,16)=mean(sparkst(ind3,7)); % tau3
                fwhm=sparkst(ind3,11); fwhm(fwhm<=0)=[];
                SPK(ssc,17)=mean(fwhm); % FWHM3
                SPK(ssc,18)=mean(sparkst(ind3,12)); % mass3
                SPK(ssc,26)=mean(sparkst(ind3,14)); % Ecc
            end
            SPK(ssc,19)=indx; SPK(ssc,20)=nspk2;
            v=[2,4,6,7,11,12];
            for ii=1:length(v)
                mm=sparkst(:,v(ii));
                ind=find(mm<=0); mm(ind)=[];
                if~isempty(mm),
                    SPKall(ssc,ii)=mean(mm);
                end
            end
            SPKall(ssc,7)=indx;
            SPKall(ssc,8)=nspk2;
            
            nRyRs(ssc,4)=size(sparkst,1);
            save([FOLc '\massspk.mat'],'massSpk','massSpk2','masst','SPK','nRyRs','nRyRs2','SPKall','sparkst');
        end
    end
end
% SPK = [amp RoR FDHM tau FWHM mass x3; drug nspk %1RyR %2RyRs % 3RyRs
%       1,7,13 ..                        19   20    21    22      23

FOLL='D:\Dades Lab\RyR-GFP_Activation_Sparks';

RF='ResRyR_12_2';                                 % name for results folder
RFo=RF;                                         % Results folder
rmw=1;                                          % Remove waves parameters
plt=0;                                          % Spark parameters plot
rs=2;                                           % ROI around ryr for measuring time signal (.4)

%%%% FREQ
SS=dir([FOLL '\*RyR-GFP*']);
factamp=6;
vsf=[1:16]; ctff=1;
for sss=1:length(vsf)%size(SS,1)
    FOLc=[FOLL '\' SS(vsf(sss)).name]; Sc=dir([FOLc '\*RyR-GFP*']);
    ind=find([Sc(:).isdir]==0); for ii=1:length(ind), Sc(ind(ii)-ii+1)=[]; end
    fwavc=zeros(size(Sc,1),4);
    Freq=zeros(size(Sc,1),8);
    
    for ssc=1:size(Sc,1) % Cell
      
            FOL=[FOLc '\' Sc(ssc).name]; S=dir([FOL '\*RyR-GFP*']);
            S([S.isdir]==0)=[];
            % for ss=1:size(S,1), imm=strcmp([S(ss).name(end-5:end-1)],'Image');
           
            freq=[]; dens=[];
             ind1=[];
            for ss=1:size(S,1),
                indls=strfind(S(ss).name,'LS');
                imgse=strfind(S(ss).name,'Image');
                sestim=strfind(S(ss).name,'Serie1');
                if any([~isempty(indls),~isempty(imgse),~isempty(sestim)])
                    ind1=[ind1;ss];
                end
            end
            if ~isempty(ind1), S(ind1)=[];end
            
            fwav=0; fspk=0; ttw=[]; facts=[]; factsa=[]; RyRsaR=[]; wav=[];
            for ss=1:size(S,1)
                folder=[FOL '\' S(ss).name];
                % drug index
                ind=find(folder=='_'); ind1=ind(end-2);ind2=ind(end-1); drug=folder(ind1+1:ind2-1); indx=0;
                if strfind(drug,'RyR')==1, indx=1;  end % CON
                if strfind(drug,'CIL')==1, indx=2;  end, if strfind(drug,'ISO')==1, indx=3; end %% CIL CIL+CPA
                if strfind(drug,'RO')==1, indx=4; end, if strfind(drug,'FENO')==1, indx=5; end % RO RO+CPA
                if strfind(drug,'RO+CIL')==1, indx=6; end, if strfind(drug,'CGS')==1, indx=7; end % RO+CIL CGS
                if indx==0, ind1=ind(end-3);ind2=ind(end-2); drug=folder(ind1+1:ind2-1); % hi ha algun comentari
                    if ~isempty(strfind(drug,'RyR')==1)||~isempty(strfind(drug,'pl2')==1), indx=1; end % CON
                    if strfind(drug,'CIL')==1, indx=2; end,  if strfind(drug,'CIL+CPA')==1, indx=3; end,  if strfind(drug,'RO')==1, indx=4; end, if strfind(drug,'RO+CPA')==1, indx=5; end, if strfind(drug,'RO+CIL')==1, indx=6; end, if strfind(drug,'CGS')==1, indx=7; end % CGS
                end
                FOLt{ctff,2}=indx; ctff=ctff+1;
                
                 % Carrega dades:
                RF=RFo; G=load([folder,'\',RF,'\finaldetection.mat']);DX=G.DX; DT=G.DT;
                RF=[RF '_Events']; tagr='ch01'; Ch=1; load([folder,'\',RF,'\volums.mat']);
                load([folder,'\',RF,'\signals.mat']);
                Q=G.Q{1}; pos=[Q(:,3),Q(:,2)]; rrr=length(pos);
                Iryr=zeros(size(G.Isum)); for ii=1:length(pos), Iryr(Q(ii,3),Q(ii,2))=1; end;
                rr=length(pos); scz = get( 0, 'Screensize' );nnr=size(V,3);
                load([folder,'\',RF,'\ev_det.mat']);
                try,load([folder,'\',RF,'\ClustEv_' num2str(factamp) '.mat']); end
                try,load([folder '\' RF '\RyRsandEv_' num2str(factamp) '.mat']);end
                try,load([folder,'\',RF,'\trajectories.mat']);end
                try, load([folder,'\',RF,'\Det3DEvents_parameters.mat']); clear L; clear L2; end
                wavi=0;
                
                fspki=[];factsi=[];RyRsactRi=[];freqactia=[];
                if ~isempty(EVfilt),try,  
                        fspki=length(unique([EVfilt{:,10}]))/(nnr*DT); 
                        factsi=size(EVfilt,1)/(nnr*DT*rr);
                        nRyRsi=rr;
                        RyRsact=length(unique([EVfilt{:,1}]));
                        RyRsactRi=RyRsact/nRyRsi;
                        Iact=sum(Vdet,3); act=Iact(:); act(act==0)=[];
                        freqactia=mean(act)/(nnr*DT*RyRsact);
                    end, end
                if ~isempty(any(eventsW>0)),
                    Ic=conv2(eventsW,ones(5,round(20*DT)),'same'); Ic(Ic>0)=1; % imagesc(Ic);
                    L=bwlabel(Ic); % imagesc(L);
                    uu=unique(L); uu(uu==0)=[];
                    
                    ll=sum(L,1);
                    tt=length(find(ll==0));
                    wavi=length(uu);
                    
                    fwavi=length(uu)/(nnr*DT*rr); % en spk/(ms*um)
                    fwav=[fwav;fwavi];
                    
                    ttw=[ttw;(nnr-tt)*DT];
                    
                   try, fspki=length(unique([EVfilt{:,10}]))/(tt*DT*rr); end % spk/(ms*um)
                   try, factsi=size(EVfilt,1)/(nnr*DT*rr); end % Freq activacions
                   try,  freqactia=mean(act)/(nnr*DT*RyRsact); end % Freq RyRs actius
                   
                end
                 fspk=[fspk;fspki];
                 facts=[facts;factsi];
                 factsa=[factsa;freqactia];
                 RyRsaR=[RyRsaR;RyRsactRi];
                 wav=[wav;wavi];
            end

        fwavc(ssc,1)=mean(fwav);
        fwavc(ssc,2)=indx;
        fwavc(ssc,3)=mean(fspk);
        if ~isempty(ttw),fwavc(ssc,4)=mean(ttw);end % save([FOLc '\freq.mat'],'fwavc');
        if any(fspk>0)||any(fwav>0),
        Freq(ssc,1)=mean(fspk);
        Freq(ssc,2)=mean(fwav);
        Freq(ssc,3)=mean(facts);
        Freq(ssc,4)=mean(factsa);
        Freq(ssc,5)=mean(RyRsaR);
        Freq(ssc,6)=indx;
        Freq(ssc,7)=sum(wav);
        Freq(ssc,8)=size(S,1)*nnr*DT*10^-3; % temps exps
        else
            Freq(ssc,6)=indx;
        end
        save([FOLc '\freq.mat'],'fwavc','Freq');
    end
    
end


clear all;   
%%%% TOTS!
FOLL='D:\Dades Lab\RyR-GFP_Activation_Sparks';
SS=dir([FOLL '\*RyR-GFP*']);
massSpkt=[]; SPKt=[]; sparksttt=[];
vsf=[1:16]; % vsf=[8:14];
ctff=1; masst=[];nRyRst=[]; Freqt=[]; FreqR=[]; massSpkt2=[];SPKallt=[];
for sss=1:size(vsf,2)
    % HISTOGRAMS:
    FOL=['D:\Dades Lab\RyR-GFP_Activation_Sparks\' SS(vsf(sss)).name];
    load([FOL '\freq.mat']);
    load([FOL '\massspk.mat']);
    massSpkt=[massSpkt;massSpk];
     massSpkt2=[massSpkt2;massSpk2];
    SPKt=[SPKt;SPK];
    SPKallt=[SPKallt;SPKall];
    nRyRst=[nRyRst;nRyRs2];
    Freqt=[Freqt;fwavc];
    FreqR=[FreqR;Freq];
end

% MRAM & MLAM (per MVM canviar la crida de carpetes)
FOLL='E:\Dades\Calcium Imaging - Calgary\RyR-GFP_MAM\MRAM';
FOLL='E:\Dades\Calcium Imaging - Calgary\RyR-GFP_MAM\MLAM';

sparksttt=[];
SS=[dir([FOLL '\*HET*'])];
factamp=6;
vsf=[1:size(SS,1)]; ctff=1; masst=[];
ctcell=1;RF='ResRyR_12_2_2';
for sss=1:length(vsf)%size(SS,1)
    FOLc=[FOLL '\' SS(vsf(sss)).name]; Sc=dir([FOLc '\*Cell*']);
    ind=find([Sc(:).isdir]==0);
    for ii=1:length(ind), Sc(ind(ii)-ii+1)=[]; end
    for ssc=1:size(Sc,1) % Cell
        if Sc(ssc).isdir,
            FOL=[FOLc '\' Sc(ssc).name]; S=[dir([FOL '\*Rest*']);dir([FOL '\*job*'])];
           
            S([S.isdir]==0)=[];
            % for ss=1:size(S,1), imm=strcmp([S(ss).name(end-5:end-1)],'Image');
            nRyRs2i=zeros(size(S,1),11);
            freq=[]; dens=[];
             ind1=[];
            for ss=1:size(S,1),
                indls=strfind(S(ss).name,'LS');
                imgse=strfind(S(ss).name,'Image');
                sestim=strfind(S(ss).name,'Serie1');
                if any([~isempty(indls),~isempty(imgse),~isempty(sestim)])
                    ind1=[ind1;ss];
                end
            end
            if ~isempty(ind1), S(ind1)=[];end
            sparkstt=[];
            for ss=1:size(S,1)
                folder=[FOL '\' S(ss).name];
                % drug index
                nameindx=FOLL(find(FOLL=='\',1,'last')+1:end);
                if strcmp(nameindx,'MRAM')==1, indx=1; else, indx=2; end;
%                 ind=find(folder=='_'); ind1=ind(end-2);ind2=ind(end-1); drug=folder(ind1+1:ind2-1); indx=0;
%                 if strfind(drug,'RyR')==1, indx=1;  end % CON
%                 if strfind(drug,'CIL')==1, indx=2;  end, if strfind(drug,'ISO')==1, indx=3; end %% CIL CIL+CPA
%                 if strfind(drug,'RO')==1, indx=4; end, if strfind(drug,'FENO')==1, indx=5; end % RO RO+CPA
%                 if strfind(drug,'RO+CIL')==1, indx=6; end, if strfind(drug,'CGS')==1, indx=7; end % RO+CIL CGS
%                 if indx==0, ind1=ind(end-3);ind2=ind(end-2); drug=folder(ind1+1:ind2-1); % hi ha algun comentari
%                     if ~isempty(strfind(drug,'RyR')==1)||~isempty(strfind(drug,'pl2')==1), indx=1; end % CON
%                     if strfind(drug,'CIL')==1, indx=2; end,  if strfind(drug,'CIL+CPA')==1, indx=3; end,  if strfind(drug,'RO')==1, indx=4; end, if strfind(drug,'RO+CPA')==1, indx=5; end, if strfind(drug,'RO+CIL')==1, indx=6; end, if strfind(drug,'CGS')==1, indx=7; end % CGS
%                 end
                ind=find(folder=='\',1,'last');
                name=folder(ind+1:end);
                 FOLt{ctff,1}=name; 
                FOLt{ctff,2}=indx; 
                FOLt{ctff,3}=ctff; ctff=ctff+1;
                load([FOLc '\Sparks_parameters.mat']);
                sparkst=[sparkst,ones(size(sparkst,1),1)*ctff];
                
            end
            sparkstt=[sparkstt;sparkst];
            
            % save([FOLc '\massspk.mat'],'massSpk','massSpk2','masst','SPK','nRyRs','nRyRs2','SPKall','sparkst');
        end
    end
     sparksttt=[sparksttt;sparkstt];
end

sparksttt(:,14)=[];
sparksttt(:,1)=[];
Vnames={'amp','AMP','RoR','t2p','FDHM','tau','r2','bl','nRyRs','FWHM','mass','Drug','number_cell'};

T=array2table(sparksttt,'VariableNames',Vnames);
writetable(T,[FOLL '\sparks_all_MLAM.xlsx'])      

%Tfol=load([FOLL '\cells_names.mat']); Tfol=Tfol.T;
Tfol=array2table(FOLt,'VariableNames',{'Cell','Drug','number_cell'});
writetable(Tfol,[FOLL '\Cells_MLAM.xlsx']) 




icon=find(massSpkt(:,3)==1);
icil=find(massSpkt(:,3)==2);
iro=find(massSpkt(:,3)==4);
iiso=find(massSpkt(:,3)==3);
ifeno=find(massSpkt(:,3)==5);


c1=[193 2 37]/255;
c2=[6 80 154]/255;
c3=[13 193 115]/255;
c4=[242 193 0]/255;
mg=.08; FS=12; LW=2;

figure(1);clf;
noms={'mass 1 RyRs','mass 2 RyRs'};
for ii=1:2
    ax=subNM(1,2,ii,mg);set(ax,'fontsize',FS);
    x1=massSpkt2(icon,ii); x2=massSpkt2(icil,ii); x3=massSpkt2(iro,ii);
    x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0); if ii==1, x1(x1>50)=0; end
    x1(isnan(x1))=[]; x2(isnan(x2))=[]; x3(isnan(x3))=[];
    CompareDistributions({x1,x3,x2},'tst','ranksum','lab',{'CON','RO','CIL'},'dsp','bar','ori',...
        'vert','tit',[noms{ii}],'sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',[c1;c4;c2]); axis square;
end
set(gcf,'position',[1 41 1920/2.5 964/3]);
saveWysiwyg(gcf,([FOLL '\mass_spk_nRyRs.png']));
% ISO FENO
figure(1);clf;
noms={'mass 1 RyRs','mass 2 RyRs'};
for ii=1:2
    ax=subNM(1,2,ii,mg);set(ax,'fontsize',FS);
    x1=massSpkt2(icon,ii); x2=massSpkt2(ifeno,ii); x3=massSpkt2(iiso,ii);
    x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0); if ii==1, x1(x1>50)=0; end
    x1(isnan(x1))=[]; x2(isnan(x2))=[]; x3(isnan(x3))=[];
    CompareDistributions({x1,x3,x2},'tst','ranksum','lab',{'CON','ISO','FENO'},'dsp','bar','ori',...
        'vert','tit',[noms{ii}],'sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',[c1;c4;c2]); axis square;
end
set(gcf,'position',[1 41 1920/2.5 964/3]);

figure(1);clf;
noms={'mass 1 RyRs','mass 2 RyRs'};
for ii=1:2
    ax=subNM(1,2,ii,mg);set(ax,'fontsize',FS);
    x1=massSpkt(icon,ii); x2=massSpkt(icil,ii); x3=massSpkt(iro,ii);
    x1(isnan(x1))=[]; x2(isnan(x2))=[]; x3(isnan(x3))=[];
     x1(x1==0)=[]; x2(x2==0)=[]; x3(x3==0)=[];
    D=[x1;x3;x2]; X=[ones(length(x1),1);2*ones(length(x3),1);3*ones(length(x2),1)];
    C=[c1;c4;c2];
    fancyBoxplot2(D,X,C,'tst','ranksum','lab',{'CON','RO','CIL'}); axis square;
    title([noms{ii}]);
    set(gca,'xtick',[1,2,3],'xticklabel',{'CON','RO','CIL'},'box','on');grid on; ylabel('Signal mass  (a.u.)'); yl=ylim();
end
set(gcf,'position',[1 41 1920/2.5 964/3]);
saveWysiwyg(gcf,([FOLL '\mass_spk_nRyRs_d.png']));

figure(1);clf;
noms={'mass 1 RyRs','mass 2 RyRs'};
% ax=subNM(1,2,1,mg);set(ax,'fontsize',FS);
x1=massSpkt(icon,1); x2=massSpkt(icil,1); x3=massSpkt(iro,1);
x12=massSpkt(icon,2); x22=massSpkt(icil,2); x32=massSpkt(iro,2);
x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
x12=x12(x12>0); x22=x22(x22>0); x32=x32(x32>0);
x1(isnan(x1))=[]; x2(isnan(x2))=[]; x3(isnan(x3))=[];
x12(isnan(x12))=[]; x22(isnan(x22))=[]; x32(isnan(x32))=[];
CompareDistributions({x1,x12,x3,x32,x2,x22},'tst','ranksum','lab',{'CON','CON 2','RO','RO 2','CIL','CIL 2'},'dsp','bar','ori',...
    'vert','tit','Mass nº RyRs spk','sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',[c1;c1*.8;c4;c4*.8;c2;c2*.8]); axis square;
set(gcf,'position',[1 41 1920/2 964/2]);
saveWysiwyg(gcf,([FOLL '\mass_spk_nRyRs_2.png']));

c11=[255, 204, 0]/255;
c55=[244,106,78]/255;
c33=[45, 89, 134]/255;

c33=[255, 204, 0]/255;
c55=[0, 134, 179]/255;
c11=[115, 0, 153]/255;


icon=find(massSpkt2(:,3)==1);
icil=find(massSpkt2(:,3)==2);
iro=find(massSpkt2(:,3)==4);
iiso=find(massSpkt2(:,3)==3);
ifeno=find(massSpkt2(:,3)==5);
% COMP SPARKS

% ISO FENO
figure(1);clf;
noms={'mass 1 RyRs','mass 2 RyRs'};
% ax=subNM(1,2,1,mg);set(ax,'fontsize',FS);
x1=massSpkt2(icon,1); x2=massSpkt2(ifeno,1); x3=massSpkt2(iiso,1);
x12=massSpkt2(icon,2); x22=massSpkt2(ifeno,2); x32=massSpkt2(iiso,2);
x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
x12=x12(x12>0); x22=x22(x22>0); x32=x32(x32>0); x32(x32<8)=[];
x1(isnan(x1))=[]; x2(isnan(x2))=[]; x3(isnan(x3))=[];
x12(isnan(x12))=[]; x22(isnan(x22))=[]; x32(isnan(x32))=[];
CompareDistributions({x1,x12,x2,x22,x3,x32},'tst','ranksum','lab',{'CON','CON 2','FENO','FENO 2','ISO','ISO 2'},'dsp','bar','ori',...
    'vert','tit','Mass nº RyRs spk','sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',[c11;c11*.8;c33;c33*.8;c55;c55*.8]); axis square;
saveWysiwyg(gcf,([FOLL '\mass_spk_nRyRs_ISO_FENO.png']));

figure(1);clf;
noms={'mass 1 RyRs','mass 2 RyRs'};
% ax=subNM(1,2,1,mg);set(ax,'fontsize',FS);
x1=massSpkt2(icon,1); x2=massSpkt2(ifeno,1); x3=massSpkt2(iiso,1);
x12=massSpkt2(icon,2); x22=massSpkt2(ifeno,2); x32=massSpkt2(iiso,2);
x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
x12=x12(x12>0); x22=x22(x22>0); x32=x32(x32>0); x32(x32<8)=[];
x1(isnan(x1))=[]; x2(isnan(x2))=[]; x3(isnan(x3))=[];
x12(isnan(x12))=[]; x22(isnan(x22))=[]; x32(isnan(x32))=[];x12(x12>60)=0;
subplot(131)
CompareDistributions({x1,x12},'tst','ranksum','lab',{'1 RyR2','2 RyR2'},'dsp','bar','ori',...
    'vert','tit','CON','sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',...
    [c11;c11*.8;],'wi',.4); axis square;
ylim([0 35]);
subplot(132)
CompareDistributions({x2,x22},'tst','ranksum','lab',{'1 RyR2','2 RyR2'},'dsp','bar','ori',...
    'vert','tit','FENO','sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',...
    [c33;c33*.8;],'wi',.4); axis square;
ylim([0 35]);
subplot(133)
CompareDistributions({x3,x32},'tst','ranksum','lab',{'1 RyR2','2 RyR2'},'dsp','bar','ori',...
    'vert','tit','ISO','sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',...
    [c55;c55*.8;],'wi',.4); axis square;
ylim([0 35]);
set(gcf,'position',[1 41 1920/3 964/3]);
saveWysiwyg(gcf,([FOLL '\mass_spk_nRyRs_ISO_FENO_2.png']));

figure(1);clf;
noms={'mass 1 RyRs','mass 2 RyRs'};
% ax=subNM(1,2,1,mg);set(ax,'fontsize',FS);
x1=massSpkt(icon,1); x2=massSpkt(ifeno,1); x3=massSpkt(iiso,1);
x12=massSpkt(icon,2); x22=massSpkt(ifeno,2); x32=massSpkt(iiso,2);
x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
x12=x12(x12>0); x22=x22(x22>0); x32=x32(x32>0); x32(x32<8)=[];
x1(isnan(x1))=[]; x2(isnan(x2))=[]; x3(isnan(x3))=[];
x12(isnan(x12))=[]; x22(isnan(x22))=[]; x32(isnan(x32))=[];% x12(x12>60)=0;
subplot(131)
CompareDistributions({x1,x2},'tst','ranksum','lab',{'CON','FENO'},'dsp','bar','ori',...
    'vert','tit','Single RyR2 spark','sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',...
    [c11;c33;],'wi',.4); axis square;
ylim([0 35]);
subplot(132)
CompareDistributions({x12,x22},'tst','ranksum','lab',{'CON','FENO'},'dsp','bar','ori',...
    'vert','tit','Double RyR2 spark','sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',...
    [c11*.8;c33*.8;],'wi',.4); axis square;
ylim([0 35]);
subplot(133)
CompareDistributions({x3,x32},'tst','ranksum','lab',{'1 RyR2','2 RyR2'},'dsp','bar','ori',...
    'vert','tit','ISO','sig','*', 'uni','Signal mass  (a.u.)','fig',0,'fs',FS,'dc',...
    [c55;c55*.8;],'wi',.4); axis square;
ylim([0 35]);
set(gcf,'position',[1 41 1920/3 964/3]);
saveWysiwyg(gcf,([FOLL '\mass_spk_nRyRs_ISO_FENO_3.png']));

% ISO FENO RyRs
figure(1);clf;
% ax=subNM(1,2,1,mg);set(ax,'fontsize',FS);
x1=nRyRst(icon,6); x2=nRyRst(ifeno,6); x3=nRyRst(iiso,6);
x12=nRyRst(icon,7); x22=nRyRst(ifeno,7); x32=nRyRst(iiso,7);
x13=nRyRst(icon,8); x23=nRyRst(ifeno,8); x33=nRyRst(iiso,8);
x14=nRyRst(icon,9); x24=nRyRst(ifeno,9); x34=nRyRst(iiso,9);
% x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
% x12=x12(x12>0); x22=x22(x22>0); x32=x32(x32>0);
% x13=x13(x13>0); x23=x23(x23>0); x33=x33(x33>0); x13(x13>0.2)=[];
% x14=x14(x14>0); x24=x24(x24>0); x34=x34(x34>0);
x1f=nRyRst(icon,11); x2f=nRyRst(ifeno,11); x3f=nRyRst(iiso,11);
x1fw=Freqt(icon,1); x2fw=Freqt(ifeno,1); x3fw=Freqt(iiso,1);
% x1f=x1f(x1f>0); x2f=x2f(x2f>0); x3f=x3f(x3f>0);
subplot(231)
CompareDistributions({x1,x2,x3},'tst','ranksum','lab',{'CON','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit','Spk 1 RyRs','sig','*', 'uni','fraction (n/N)','fig',0,...
    'fs',FS,'dc',[c11;c33;c55]); axis square;
subplot(232)
CompareDistributions({x12,x22,x32},'tst','ranksum','lab',{'CON','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit','Spk 2 RyRs','sig','*', 'uni','fraction (n/N)','fig',0,...
    'fs',FS,'dc',[c11;c33;c55]); axis square;
subplot(233)
CompareDistributions({x13,x23,x33},'tst','ranksum','lab',{'CON','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit','Spk 3 RyRs','sig','*', 'uni','fraction (n/N)','fig',0,...
    'fs',FS,'dc',[c11;c33;c55]); axis square;
subplot(234)
% CompareDistributions({x14,x24,x34},'tst','ranksum','lab',{'CON','FENO','ISO'},'dsp','bar','ori',...
%     'vert','tit','Spk >3 RyRs','sig','*', 'uni','fraction (n/N)','fig',0,...
%     'fs',FS,'dc',[c11;c33;c55]); axis square;
CompareDistributions({x1f*100,x2f*100,x3f*100},'tst','ranksum','lab',{'CON','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit','Spk frequency','sig','*', 'uni','sparks/(s*\mum)','fig',0,...
    'fs',FS,'dc',[c11;c33;c55]); axis square;
subplot(235)
CompareDistributions({x1fw*100,x2fw*100,x3fw*100},'tst','ranksum','lab',{'CON','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit','Waves frequency','sig','*', 'uni','waves/s','fig',0,...
    'fs',FS,'dc',[c11;c33;c55]); axis square;
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOLL '\spk_nRyRs_ISO_FENO.png']));

% ISO FENO RyRs
ind1=find(SPKt(:,22)==0); % Troba cells amb spks només d'un RyR
ind2=find(SPKt(:,22)>0); % cells amb sparks de 2 RyRs
cc=intersect(icon,ind1);
cc2=intersect(icon,ind2);
cc3=intersect(ifeno,ind1);
cc4=intersect(ifeno,ind2);
figure(1);clf;
% ax=subNM(1,2,1,mg);set(ax,'fontsize',FS);
x1=nRyRst(icon,6); x2=nRyRst(ifeno,6); x3=nRyRst(iiso,6);
x12=nRyRst(icon,7); x22=nRyRst(ifeno,7); x32=nRyRst(iiso,7);
x13=nRyRst(icon,8); x23=nRyRst(ifeno,8); x33=nRyRst(iiso,8);
x14=nRyRst(icon,9); x24=nRyRst(ifeno,9); x34=nRyRst(iiso,9);
% x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
% x12=x12(x12>0); x22=x22(x22>0); x32=x32(x32>0);
% x13=x13(x13>0); x23=x23(x23>0); x33=x33(x33>0); x13(x13>0.2)=[];
% x14=x14(x14>0); x24=x24(x24>0); x34=x34(x34>0);
x1f=nRyRst(icon,11); x2f=nRyRst(ifeno,11); x3f=nRyRst(iiso,11);
x1fw=Freqt(icon,1); x2fw=Freqt(ifeno,1); x3fw=Freqt(iiso,1);
x11=FreqR(cc,1); x122=FreqR(cc2,1); 
x1f1=FreqR(cc3,1); x1f2=FreqR(cc4,1); 
% x1f=x1f(x1f>0); x2f=x2f(x2f>0); x3f=x3f(x3f>0);
subplot(231)
CompareDistributions({x1,x2},'tst','ranksum','lab',{'CON','FENO'},'dsp','bar','ori',...
    'vert','tit','Sparks 1 RyR2','sig','*', 'uni','fraction','fig',0,...
    'fs',FS,'dc',[c11;c33],'wi',0.4); axis square;
subplot(232)
CompareDistributions({x12,x22},'tst','ranksum','lab',{'CON','FENO'},'dsp','bar','ori',...
    'vert','tit','Sparks 2 RyR2','sig','*', 'uni','fraction','fig',0,...
    'fs',FS,'dc',[c11;c33],'wi',0.4); axis square;
subplot(233)
CompareDistributions({x13,x23},'tst','ranksum','lab',{'CON','FENO'},'dsp','bar','ori',...
    'vert','tit','Sparks 3 RyR2','sig','*', 'uni','fraction','fig',0,...
    'fs',FS,'dc',[c11;c33],'wi',0.4); axis square;
subplot(234)
% CompareDistributions({x14,x24,x34},'tst','ranksum','lab',{'CON','FENO','ISO'},'dsp','bar','ori',...
%     'vert','tit','Spk >3 RyRs','sig','*', 'uni','fraction (n/N)','fig',0,...
%     'fs',FS,'dc',[c11;c33;c55]); axis square;
CompareDistributions({x1f*100,x2f*100},'tst','ranksum','lab',{'CON','FENO'},'dsp','bar','ori',...
    'vert','tit','Spk frequency','sig','*', 'uni','sparks/(s·\mum)','fig',0,...
    'fs',FS,'dc',[c11;c33],'wi',0.4); axis square;
subplot(235)
CompareDistributions({x1fw*100,x2fw*100},'tst','ranksum','lab',{'CON','FENO'},'dsp','bar','ori',...
    'vert','tit','Waves frequency','sig','*', 'uni','waves/(s·\mum)','fig',0,...
    'fs',FS,'dc',[c11;c33],'wi',0.4); axis square;
set(gcf,'position',[1 41 1920/2.5 964/2]);
subplot(236)
CompareDistributions({x11,x1f1,x122,x1f2},'tst','ranksum','lab',{'C1','F1','C2','F2'},'dsp','bar','ori',...
    'vert','tit','sparks frequency','sig','*', 'uni','waves/s','fig',0,...
    'fs',FS,'dc',[c11;c33;c11;c33;]); axis square;
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOLL '\spk_nRyRs_FENO.png']));

std(x1)/sqrt(length(x1))
std(x2)/sqrt(length(x2))

std(x12)/sqrt(length(x12))
std(x22)/sqrt(length(x22))

std(x13)/sqrt(length(x13))
std(x23)/sqrt(length(x23))

std(x1f)/sqrt(length(x1f))
std(x2f)/sqrt(length(x2f))

std(x1fw)/sqrt(length(x1fw))
std(x2fw)/sqrt(length(x2fw))



figure(1);clf;
noms={'mass 1 RyRs','mass 2 RyRs'};
% ax=subNM(1,2,1,mg);set(ax,'fontsize',FS);
x1=massSpkt(icon,1); x2=massSpkt(icil,1); x3=massSpkt(iro,1);
x12=massSpkt(icon,2); x22=massSpkt(icil,2); x32=massSpkt(iro,2);
x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
x12=x12(x12>0); x22=x22(x22>0); x32=x32(x32>0);
x1(isnan(x1))=[]; x2(isnan(x2))=[]; x3(isnan(x3))=[];
x12(isnan(x12))=[]; x22(isnan(x22))=[]; x32(isnan(x32))=[];
C=[c1;c1*.8;c4;c4*.8;c2;c2*.8];
D=[x1;x12;x3;x32;x2;x22];
X=[ones(length(x1),1);2*ones(length(x12),1);3*ones(length(x3),1);...
    4*ones(length(x32),1);5*ones(length(x2),1);6*ones(length(x22),1)];
fancyBoxplot2(D,X,C,'tst','ranksum','lab',{'CON','CON 2','RO','RO 2','CIL','CIL 2'}); axis square;
title('nº RyRs');
set(gca,'xtick',[1,2,3,4,5,6],'xticklabel',{'CON','CON 2','RO','RO 2','CIL','CIL 2'},'box','on');
grid on; ylabel('Signal mass  (a.u.)'); yl=ylim();
set(gcf,'position',[1 41 1920/2 964/2]);
saveWysiwyg(gcf,([FOLL '\mass_spk_nRyRs_2_distrb.png']));

% Freq
icon=find(FreqR(:,6)==1);
icil=find(FreqR(:,6)==2);
iro=find(FreqR(:,6)==4);
ifeno=find(FreqR(:,6)==5);
iiso=find(FreqR(:,6)==3);

figure(1);clf;
v=[3 4 5 1 2];
noms={'Activations freq','Activ. Active freq','RyRs act. Ratio','Spk freq','Waves freq'};
uni={'sparks/(s·ROI)','sparks/(s·ROI actives)','#RyRs act/#RyRs','sparks/(s·\mum)','waves/(s·\mum)'};
for ii=1:5
    x1=FreqR(icon,v(ii));
    x2=FreqR(ifeno,v(ii));
    x3=FreqR(iiso,v(ii));
    subplot(2,3,ii)
    CompareDistributions({x1,x2,x3},'tst','ranksum','lab',{'CON','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{ii},'sig','*', 'uni',uni{ii},'fig',0,...
    'fs',FS,'dc',[c11;c33;c55]); axis square;
end
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOLL '\spk_nRyRs_ISO_FENO.png']));





% Sparks morfology
[ind,~]=find(isnan(SPKt)); SPKt(ind,:)=[]; % Fora els NaN

icon=find(SPKt(:,19)==1);
icil=find(SPKt(:,19)==2);
iro=find(SPKt(:,19)==4);
ifeno=find(SPKt(:,19)==5);
iiso=find(SPKt(:,19)==3);
indd={icon,icil,iro}; cc=[c1;c2;c4];
namesd={'CON','CIL','RO'};
noms={'amp','RoR','FDHM','tau','FWHM','mass'};
xl={'F/Fo','1/ms','ms','ms','um','(a.u.)'};
figure(1);
for ii=1:3 % drugs
    ind=indd{ii}; clf;
    for jj=1:6
        ax=subNM(2,3,jj,mg);set(ax,'fontsize',FS);
        
        x1=SPKt(ind,jj); x2=SPKt(ind,jj+6);x3=SPKt(ind,jj+12);
        x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
        CompareDistributions({x1,x2,x3},'tst','ranksum','lab',{'1','2','3'},'dsp','bar','ori',...
            'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,'dc',[cc(ii,:);cc(ii,:);cc(ii,:)]);
        axis square;
    end
    % title([namesd{ii}]);
    set(gcf,'position',[1 41 1920/1.8 964/1.8]);
    saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_' namesd{ii} '.png']));
end

icon=find(SPKt(:,19)==1);
icil=find(SPKt(:,19)==2);
iro=find(SPKt(:,19)==4);
ifeno=find(SPKt(:,19)==5);
iiso=find(SPKt(:,19)==3);
indd={icon,ifeno,iiso}; cc=[c11;c33;c55];
namesd={'CON','FENO','ISO'};
noms={'amp','RoR','FDHM','tau','FWHM','mass'};
xl={'F/Fo','1/ms','ms','ms','um','(a.u.)'};
figure(1);
for ii=1:3 % drugs
    ind=indd{ii}; clf;
    for jj=1:6
        ax=subNM(2,3,jj,mg);set(ax,'fontsize',FS);
        
        x1=SPKt(ind,jj); x2=SPKt(ind,jj+6);x3=SPKt(ind,jj+12);
        x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
        CompareDistributions({x1,x2,x3},'tst','ranksum','lab',{'1','2','3'},'dsp','bar','ori',...
            'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,'dc',[cc(ii,:);cc(ii,:);cc(ii,:)]);
        axis square;
    end
    % title([namesd{ii}]);
    set(gcf,'position',[1 41 1920/1.8 964/1.8]);
    saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_' namesd{ii} '_2.png']));
end

namesd={'CON','FENO','ISO'};
noms={'amp','RoR','FDHM','tau','FWHM','mass'};
xl={'F/Fo','1/ms','ms','ms','um','(a.u.)'};
figure(1); clf;
% cell 43 CON!
for jj=1:3%6
    ax=subNM(1,3,jj,mg);set(ax,'fontsize',FS);
    
    x1=SPKt(icon,jj); x2=SPKt(icon,jj+6);x3=SPKt(icon,jj+12);
    
    x1f=SPKt(ifeno,jj); x2f=SPKt(ifeno,jj+6);x3f=SPKt(ifeno,jj+12);
    x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);% x2(8)=[];
   
    x1f=x1f(x1f>0); x2f=x2f(x2f>0); x3f=x3f(x3f>0);
    CompareDistributions({x1,x2,x1f,x2f},'tst','ranksum','lab',{'C1','C2','F1','F2'},'dsp','bar','ori',...
        'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,...
        'dc',[cc(1,:);cc(1,:);cc(2,:);cc(2,:);],'wi',0.4);
    axis square;
end
% title([namesd{ii}]);
set(gcf,'position',[1 41 1920/3 964/3]);
saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_FENO_1.png']));

namesd={'CON','FENO','ISO'};
noms={'amp','RoR','FDHM','tau','FWHM','mass'};
xl={'F/Fo','1/ms','ms','ms','um','(a.u.)'};
figure(1); clf;
% cell 43 CON!
for jj=1:3%6
    ax=subNM(1,3,jj,mg);set(ax,'fontsize',FS);
    
    x1=SPKt(icon,jj); x2=SPKt(icon,jj+6);x3=SPKt(icon,jj+12);
    
    x1f=SPKt(ifeno,jj); x2f=SPKt(ifeno,jj+6);x3f=SPKt(ifeno,jj+12);
    x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);% x2(8)=[];
   
    x1f=x1f(x1f>0); x2f=x2f(x2f>0); x3f=x3f(x3f>0);
    CompareDistributions({x1,x1f,[],x2,x2f},'tst','ranksum','lab',{'C1','F1','','C2','F2'},'dsp','bar','ori',...
        'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,...
        'dc',[cc(1,:);cc(2,:);cc(1,:);cc(1,:);cc(2,:);],'wi',0.5);
    axis square;
end
% title([namesd{ii}]);
set(gcf,'position',[1 41 1920/3 964/3]);
saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_FENO_12.png']));

mean(x1)
std(x1)/sqrt(length(x1))
mean(x1f)
std/sqrt(length(x1f))
mean(x2)/sqrt(std(x2))
mean(x2f)/sqrt(std(x2f))


namesd={'CON','FENO','ISO'};
noms={'amp','RoR','FDHM','tau','FWHM','mass'};
xl={'F/Fo','1/ms','ms','ms','um','(a.u.)'};
figure(1); clf;
% cell 43 CON!
for jj=4:5%6
     ax=subNM(1,3,jj-2,mg);set(ax,'fontsize',FS);
    
    x1=SPKt(icon,jj); x2=SPKt(icon,jj+6);x3=SPKt(icon,jj+12);
    
    x1f=SPKt(ifeno,jj); x2f=SPKt(ifeno,jj+6);x3f=SPKt(ifeno,jj+12);
    x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0); %x2(8)=[];
   
    x1f=x1f(x1f>0); x2f=x2f(x2f>0); x3f=x3f(x3f>0);
    CompareDistributions({x1,x2,x1f,x2f},'tst','ranksum','lab',{'C1','C2','F1','F2'},'dsp','bar','ori',...
        'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,...
        'dc',[cc(1,:);cc(1,:);cc(2,:);cc(2,:);],'wi',0.4);
    axis square;
end
% title([namesd{ii}]);
set(gcf,'position',[1 41 1920/3 964/3]);
saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_FENO_2.png']));


namesd={'CON','FENO','ISO'};
noms={'amp','RoR','FDHM','tau','FWHM','mass'};
xl={'F/Fo','1/ms','ms','ms','um','(a.u.)'};
figure(1); clf;
% cell 43 CON!
for jj=4:5%6
     ax=subNM(1,3,jj-2,mg);set(ax,'fontsize',FS);
    
    x1=SPKt(icon,jj); x2=SPKt(icon,jj+6);x3=SPKt(icon,jj+12);
    
    x1f=SPKt(ifeno,jj); x2f=SPKt(ifeno,jj+6);x3f=SPKt(ifeno,jj+12);
    x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0); %x2(8)=[];
   
    x1f=x1f(x1f>0); x2f=x2f(x2f>0); x3f=x3f(x3f>0);
      CompareDistributions({x1,x1f,[],x2,x2f},'tst','ranksum','lab',{'C1','F1','','C2','F2'},'dsp','bar','ori',...
        'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,...
        'dc',[cc(1,:);cc(2,:);cc(1,:);cc(1,:);cc(2,:);],'wi',0.5);
    axis square;
end
% title([namesd{ii}]);
set(gcf,'position',[1 41 1920/3 964/3]);
saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_FENO_22.png']));

namesd={'CON','FENO','ISO'};
noms={'amp','RoR','FDHM','tau','FWHM','mass'};
xl={'F/Fo','1/ms','ms','ms','um','(a.u.)'};
figure(1); clf;
% cell 43 CON!
for jj=1:6
    ax=subNM(2,3,jj,mg);set(ax,'fontsize',FS);
    
    x1=SPKt(icon,jj); x2=SPKt(icon,jj+6);x3=SPKt(icon,jj+12);
    x1i=SPKt(iiso,jj); x2i=SPKt(iiso,jj+6);x3i=SPKt(iiso,jj+12);
    x1f=SPKt(ifeno,jj); x2f=SPKt(ifeno,jj+6);x3f=SPKt(ifeno,jj+12);
    x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
    x1i=x1i(x1i>0); x2i=x2i(x2i>0); x3i=x3i(x3i>0);
    x1f=x1f(x1f>0); x2f=x2f(x2f>0); x3f=x3f(x3f>0);
    
%     x3(2)=[]; x3(5)=[];
%     x3i(1)=[]; x3i(4)=[];
%     x3f(5)=[];
    CompareDistributions({x1,x2,x3,x1i,x2i,x3i,x1f,x2f,x3f},'tst','ranksum','lab',...
        {'C1','C2','C3','I1','I2','I3','F1','F2','F3'},'dsp','bar','ori',...
        'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,'dc',...
        [cc(1,:);cc(1,:);cc(1,:);cc(3,:);cc(3,:);cc(3,:);cc(2,:);cc(2,:);cc(2,:);]);
    axis square;
end
% title([namesd{ii}]);
set(gcf,'position',[1 41 1920/1.8 964/1.8]);
saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_ISO_FENO_3.png']));




% CON nRyRs
noms={'amp','RoR','FDHM','tau','FWHM','mass'};
xl={'F/Fo','1/ms','ms','ms','um','(a.u.)'};
cc2=[6 126 186; 228 144 0]/255;
cc2=[.6 .6 .6; 0 0 0];
parametres=[1,3,5,6];
figure(1);
for ii=1:1 % drugs
    ind=indd{ii}; clf;
    for jjp=1:length(parametres)
        jj=parametres(jjp);
        ax=subNM(1,4,jjp,mg);set(ax,'fontsize',FS);
        
        x1=SPKt(ind,jj); x2=SPKt(ind,jj+6);%x3=SPKt(ind,jj+12);
        x1=x1(x1>0); x2=x2(x2>0); %x3=x3(x3>0);
        CompareDistributions({x1,x2},'tst','ranksum','lab',{'1','2'},'dsp','bar','ori',...
            'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,'dc',[cc2(1,:);cc2(2,:);],'ec',[0 0 0;0 0 0]);
        axis square;
    end
    % title([namesd{ii}]);
    set(gcf,'position',[1 41 2000/1.5 964/2]);
    saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_CON_2.png']));
end
figure(1);
cc2=[6 126 186; 228 144 0; 186 42 6]/255;
cc2=[.8 .8 .8;.4 .4 .4; 0 0 0];
for ii=1:1 % drugs
    ind=indd{ii}; clf;
    for jjp=1:length(parametres)
        jj=parametres(jjp);
        ax=subNM(1,4,jjp,mg);set(ax,'fontsize',FS);
        
        x1=SPKt(ind,jj); x2=SPKt(ind,jj+6);x3=SPKt(ind,jj+12);
        x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
        CompareDistributions({x1,x2,x3},'tst','ranksum','lab',{'1','2','3'},'dsp','bar','ori',...
            'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,'dc',[cc2(1,:);cc2(2,:);cc2(3,:);],'ec',[0 0 0;0 0 0;0 0 0]);
        axis square;
    end
    % title([namesd{ii}]);
    set(gcf,'position',[1 41 2000 964/2]);
    saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_CON_3.png']));
end

figure(1);
for ii=2:3 % drugs
    ind=indd{ii}; clf;
    for jj=1:6
        ax=subNM(2,3,jj,mg);set(ax,'fontsize',FS);
        x1c=SPKt(icon,jj); x2c=SPKt(icon,jj+6);x3c=SPKt(icon,jj+12);
        x1=SPKt(ind,jj); x2=SPKt(ind,jj+6);x3=SPKt(ind,jj+12);
        x1c=x1c(x1c>0); x2c=x2c(x2c>0); x3c=x3c(x3c>0);
        x1=x1(x1>0); x2=x2(x2>0); x3=x3(x3>0);
        CompareDistributions({x1c,x2c,x3c,x1,x2,x3},'tst','ranksum','lab',{'1','2','3','1','2','3'},'dsp','bar','ori',...
            'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,'dc',[c1;c1;c1;cc(ii,:);cc(ii,:);cc(ii,:)]);
        axis square;
    end
    % title([namesd{ii}]);
    set(gcf,'position',[1 41 1920/1.5 964/1.5]);
    saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_CON_' namesd{ii} '.png']));
end

% 1 i 2 només!
figure(1);
for ii=2:3 % drugs
    ind=indd{ii}; clf;
    for jj=1:6
        ax=subNM(2,3,jj,mg);set(ax,'fontsize',FS);
        x1c=SPKt(icon,jj); x2c=SPKt(icon,jj+6);x3c=SPKt(icon,jj+12);
        x1=SPKt(ind,jj); x2=SPKt(ind,jj+6);%x3=SPKt(ind,jj+10); x2=[x2;x3];
        x1c=x1c(x1c>0); x2c=x2c(x2c>0); %x3c=x3c(x3c>0);  x2c=[x2c;x3c];
        x1=x1(x1>0); x2=x2(x2>0);% x3=x3(x3>0);
        CompareDistributions({x1c,x2c,x1,x2},'tst','ranksum','lab',{'1','>1','1','>1'},'dsp','bar','ori',...
            'vert','tit',[noms{jj}],'sig','*', 'uni',xl{jj},'fig',0,'fs',FS,'dc',[c1;c1;cc(ii,:);cc(ii,:)]);
        axis square;
    end
    % title([namesd{ii}]);
    set(gcf,'position',[1 41 1920/1.5 964/1.1]);
    saveWysiwyg(gcf,([FOLL '\Sparks_morfology_nRyRs_CON_' namesd{ii} '.png']));
end



% Histogrames com HAM
counts=[0:1:max(masst(:,1))];
% masst(masst>30)=[];
[h,c]=histc(masst(:,1),counts);
plot(counts,h/length(masst));

% % Poisson mix model
% fit(masst(:,1),'Poisson');
% ind=find(masst(:,1)<=0); masst(ind,:)=[];
% [pdca,gn] = fitdist(masst(:,1),'Normal','By',masst(:,2));
% [pdca] = fitdist(masst(:,1),'Poisson','By',masst(:,2));
% [pdca] = fitdist(masst(:,1),'Poisson');
% [lambdahat,lambdaci] = poissfit(masst(:,1));

% GM model comprovar dues poblacions:
% 1 gaussiana
gm = fitgmdist(masst(:,1),1);
AIC1= gm.AIC;
% 2 gaussianes
gm = fitgmdist(masst(:,1),2);
AIC2= gm.AIC;
% 3 gaussianes
gm = fitgmdist(masst(:,1),3);
AIC3= gm.AIC;
% % 4 gaussianes
% gm = fitgmdist(masst,4);
% AIC4= gm.AIC;
% AIC2 ha de ser més petit que AIC1, vol dir que el model és millor!

% Representa
figure(1);clf;
x = counts; x=x';
clf;
color=[255 153 51;51 153 5;51 153 255;51 2 255]./255; hold on;
plot(counts,h/length(masst),'k','LineWidth',2);
for jj = 1:gm.NumComponents
    norm = gm.ComponentProportion(jj)*normpdf(x,gm.mu(jj,1),sqrt(gm.Sigma(jj)));
    plot(x,norm,'Color',color(jj,:),'LineWidth',1)
    xx1 = gm.mu(jj); y1 = max(norm)+.01;
    txt = strcat('mu = ',num2str(gm.mu(jj)),' / Proportion = ',num2str(gm.ComponentProportion(jj)));
    text(xx1,y1,txt,'FontSize',10,'Color',color(jj,:));
end
ylabel('Fraction'); xlabel('Mass');
saveWysiwyg(gcf,([FOLL '\Mass_GMM_MVM.png']));


% 1 i 2
ind1=find(masst(:,2)==1);
ind2=find(masst(:,2)>1);

cc2=[6 126 186; 228 144 0]/255;
c11=[6 126 186]/255;
c22=[228 144 0]/255;

counts=[0:1.5:max(masst(:,1))];
figure(1);clf;
% masst(masst>30)=[];
[h1,c111]=histc(masst(ind1,1),counts);
[h2,c222]=histc(masst(ind2,1),counts);
plot(counts,h1/length(ind1),'LineWidth',2,'Color',c11); hold on;
plot(counts,h2/length(ind2),'LineWidth',2,'Color',c22);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 i 2 CON
ind1=find(masst(:,2)==1);
ind2=find(masst(:,2)>1);
indcon=find(masst(:,3)==1);

ind11=intersect(ind1,indcon);
ind22=intersect(ind2,indcon);

c11=[102 0 0]/255;
c22=[255 51 51]/255;

indices = crossvalind('HoldOut',length(indcon),.3);
train = find(indices==1);
test = find(indices==0);

counts=[0:4:max(masst(:,1))];
% masst(masst>30)=[];
[h1,c111]=histc(masst(ind11,1),counts);
[h2,c222]=histc(masst(ind22,1),counts);
[p,h,stats]=ranksum(masst(ind11,1),masst(ind22,1));
figure(1), clf;
plot(counts,h1/length(ind11),'LineWidth',2,'Color',c11); hold on;
plot(counts,h2/length(ind22),'LineWidth',2,'Color',c22);
plot(0,'k');
ylabel('fraction (counts/n)');xlabel('signal mass (a.u)');
legend({'1 RyR cluster','2 RyR clusters',['ranksum p-value ' num2str(p)]});
title('Signal mass MVM CON');
xlim([0 100]);
saveWysiwyg(gcf,([FOLL '\Mass_hist_CON_1_2.png']));

m1=median(masst(ind11,1));
m2=median(masst(ind22,1));

% GMM MODEL
counts=[0:1:max(masst(:,1))];
% masst(masst>30)=[];
[h,c]=histc(masst(indcon(train),1),counts);
% GM model comprovar dues poblacions:
% 1 gaussiana
gm = fitgmdist(masst(indcon,1),1);
AIC1= gm.AIC;
% 2 gaussianes
gm = fitgmdist(masst(indcon(train),1),2);
AIC2= gm.AIC;
% 3 gaussianes
gm = fitgmdist(masst(indcon(train),1),3);
AIC3= gm.AIC;

% Representa
x = counts; x=x';
clf;
color=[255 153 51;51 153 5;51 153 255;51 2 255]./255; hold on;
plot(counts,h/length(indcon),'k','LineWidth',2);
for jj = 1:gm.NumComponents
    norm = gm.ComponentProportion(jj)*normpdf(x,gm.mu(jj,1),sqrt(gm.Sigma(jj)));
    plot(x,norm,'Color',color(jj,:),'LineWidth',1)
    xx1 = gm.mu(jj); y1 = max(norm)+.01;
    txt = strcat('mu = ',num2str(gm.mu(jj)),' / Proportion = ',num2str(gm.ComponentProportion(jj)));
    text(xx1,y1,txt,'FontSize',10,'Color',color(jj,:));
end
ylabel('Fraction'); xlabel('Mass (a.u.)'); title('Mass Signal CON MVM');
saveWysiwyg(gcf,([FOLL '\Mass_GMM_MVM_CON.png']));

masscon=masst(indcon(test),1); lab=masst(indcon(test),2); lab(lab>1)=2;
idx = cluster(gm,masscon); idx=2-idx; idx(idx>0)=2; idx(idx==0)=1;
masscon(:,2)=lab;
masscon(:,3)=idx;


indx1=find(lab==1);
indx2=find(lab==2);

dd=idx-lab;
% 1 RyRs
ind11=find(dd(indx1)==0); % TP
ind1n=find(dd(indx2)==-1); % FP
ind1=find(dd(indx1)==1); % FN
indtn1=length(lab)-(length(ind11)+length(ind1n)+length(ind1)); % TN
TP=length(ind11); FN = length(ind1);
FP=length(ind1n); TN=indtn1;
% Predict 1 RyR:
TPR1=TP/(TP+FN); % TPR = TP / (TP+FN)
FPR1=FP/(FP+TN); % FPR = FP / (FP + TN)
ACC1=(TP+TN)/(TP+TN+FP+FN); % ACC = (TP + TN) / (TP + TN + FP + FN)

% 2 RyRs
ind11=find(dd(indx2)==0); % TP
ind1n=find(dd(indx2)==-1); % FN
ind1=find(dd(indx1)==1); % FP
indtn1=length(lab)-(length(ind11)+length(ind1n)+length(ind1)); % TN
TP=length(ind11); FP = length(ind1);
FN=length(ind1n); TN=indtn1;
% Predict 2 RyR:
TPR2=TP/(TP+FN); % TPR = TP / (TP+FN)
FPR2=FP/(FP+TN); % FPR = FP / (FP + TN)
ACC2=(TP+TN)/(TP+TN+FP+FN); % ACC = (TP + TN) / (TP + TN + FP + FN)



%%%%%%% 1 i 2 CIL
ind1=find(masst(:,2)==1);
ind2=find(masst(:,2)>1);
indcil=find(masst(:,3)==2);

ind11=intersect(ind1,indcil);
ind22=intersect(ind2,indcil);

c11=[3 39 110]/255;
c22=[97 134 206]/255;

counts=[0:4:max(masst(:,1))];
% masst(masst>30)=[];
[h1,c111]=histc(masst(ind11,1),counts);
[h2,c222]=histc(masst(ind22,1),counts);
[p,h,stats]=ranksum(masst(ind11,1),masst(ind22,1));
figure(1), clf;
plot(counts,h1/length(ind11),'LineWidth',2,'Color',c11); hold on;
plot(counts,h2/length(ind22),'LineWidth',2,'Color',c22);
plot(0,'k');
ylabel('fraction (counts/n)');xlabel('signal mass (a.u)');
legend({'1 RyR cluster','2 RyR clusters',['ranksum p-value ' num2str(p)]});
title('Signal mass MVM CIL');
xlim([0 100]);
saveWysiwyg(gcf,([FOLL '\Mass_hist_CIL_1_2.png']));

%%%%%%% 1 i 2 RO
ind1=find(masst(:,2)==1);
ind2=find(masst(:,2)>1);
indro=find(masst(:,3)==4);

ind11=intersect(ind1,indro);
ind22=intersect(ind2,indro);

c11=[139 110 14]/255;
c22=[255 210 65]/255;

counts=[0:4:max(masst(:,1))];
% masst(masst>30)=[];
[h1,c111]=histc(masst(ind11,1),counts);
[h2,c222]=histc(masst(ind22,1),counts);
[p,h,stats]=ranksum(masst(ind11,1),masst(ind22,1));
figure(1), clf;
plot(counts,h1/length(ind11),'LineWidth',2,'Color',c11); hold on;
plot(counts,h2/length(ind22),'LineWidth',2,'Color',c22);
plot(0,'k');
ylabel('fraction (counts/n)');xlabel('signal mass (a.u)');
legend({'1 RyR cluster','2 RyR clusters',['ranksum p-value ' num2str(p)]});
title('Signal mass MVM RO');
xlim([0 100]);
saveWysiwyg(gcf,([FOLL '\Mass_hist_RO_1_2.png']));


