close all;
clear all;
clc;

%%%%%
% RyRactivation_stpau.m
% mass_sparks.m
% num_dx_totals.m
% trajectories_all.m

% ACTIVATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Activation Simple %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usuari=1;           % Program requires user interaction
% General settings
N=1;                % Number of experiments to process (set N=0 to analyze all sub-folders in folder)
tagRyR='ch00';      % relevant channel identifier
RF='ResRyR_12_2';   % name for results folder % 'ResRyRt_100_2' (ResRyR_n015)
sk=0;               % skip detection
zoom=0;             % Zoom image

% RyR2 Detection settings
th=.01;             % noise factor
NI=1;               % Number of images to process
fr=450;             % number of frames to join (56)
tn=1;               % total number of resulting images
mov=0;              % detect RyRs along the time
pg=2.5;             % Pealing rings (um)

% Filtering RyR parameters
roiR=.5;            % roi radius around RyRs candidate center (in um)
thI=.25;            % intensity threshold % .25
minR=.06;           % minimum radius (in um) % .06
maxR=0.8;           % maximum radius (in um) % 0.8
sr=1;               % show rejected sparks in output
rp=0;               % Figures RyRs signals

% Event Detection
rs=.25;              % ROI around ryr for measuring time signal (.4)
tww=100;%150;              % ms 200
Tww=1000;              % ms
umbw=0.15;             % umbral non-propagationclose all %8
umbwW=0.8;           % umbral waves 10 % 6 %30
% limamp=0;

% Event grouping
wt=50;              % Temporal window (ms) wt=20; 181004 wt = 50;
ccorr=0.7;          % Correlation
dspk=2.5;%1;           % um

% Filtering parametes
zzRyR=0.33;               % Accepted RyRs mask (no RyRs in the borders accepted) in um
ampF=.36;               % F/Fo .6 1.2 181120 .1
factampt=6;      % Fact multiply bl
tauF=[6 200];          % ms 5 100
rorF=[0.005  1.5];     % 1/ms
t2pF=[3 100];          % ms 3 100
dhmF=[12 250];            % ms (6 190307)
r2F=.05;                % .29
plt=1;                  % Graphs sparks parameters
rmw=1;                  % Remove waves

% Reprocess
ppryr=0;            % RyRs detection
vid=0;              % Ca channel videos
ppsign=1;           % Signals
ppdet=1;          % Event paramerets
ppefcl=1;             % Filtering events and clustering
ppVol=1;            % Volume parameters
ppVdet=0;           % Detection 3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if usuari==1,
    FOL2=uigetdir([],'Select the main experiments folder');
else
    FOL2='D:\Dades Lab\RyR-GFP_Activation_Sparks'; mam=0;
end
Sf=dir([FOL2 '\*RyR-GFP*']);

RFo=RF;
v=[1:size(Sf,1)];

umbwWo=umbwW; 
factamp=factampt(1);
for ssf=1:length(v) % Mouse
    if usuari==1,
        FOLc=uigetdir([FOL2],'Select the main experiments folder'); Sc=dir([FOLc '\*RyR-GFP*']);
    else
        FOLc=[FOL2 '\' Sf(v(ssf)).name]; Sc=dir([FOLc '\*RyR-GFP*']);
    end
    % FOLc=[FOL2 '\' Sf(v(ssf)).name]; Sc=dir([FOLc '\*RyR-GFP*']);
    for ssc=1:size(Sc,1) % Cell
        if usuari==1,
            FOL=uigetdir([FOLc],'Select the main experiments folder'); S=dir([FOL '\*RyR-GFP*']);
        else
            FOL=[FOLc '\' Sc(ssc).name]; S=dir([FOL '\*RyR-GFP*']);
        end
        % FOL=[FOLc '\' Sc(ssc).name]; S=dir([FOL '\*RyR-GFP*']);
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
        
        for ss=1:size(S,1) % Serie
            if usuari==1,
                FOL=uigetdir([FOLc],'Select the main experiments folder'); S=dir([FOL '\*RyR-GFP*']);
            else
                FOL=[FOLc '\' Sc(ssc).name]; % S=dir([FOL '\*RyR-GFP*']);
            end
            
            folder=[FOL '\' S(ss).name];%folder=[S{ss}];
            RF=RFo;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RyRs %%%%%%%%%%%%%%%%%%%%%%%%%%
            if ppryr==1,
                ensenya(['Detecting RyRs.']); thIo=thI; fro=fr;
                if v(ssf)==12, if ssc>5, thI=0.5; fr=1000; end, end
                RyRSimple(N,tagRyR,RF,sk,th,fr,tn,pg,roiR,thI,minR,maxR,sr,folder,[],zoom,NI);
                G=load([folder,'\',RF,'\finaldetection.mat']);
                DX=G.DX; DT=G.DT;   thI=thIo; fr=fro;% [DX,DT,DZ] = buscaResolucio(folder,tagRyR)
                if mov==1, [mov,tn]=movRyRsvideo(DX, folder,tn); end
            else
                ensenya(['Loading RyRs.']);
                G=load([folder,'\',RF,'\finaldetection.mat']);
                DX=G.DX; DT=G.DT;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VID, SIGNALS & DETECTION EVENTS
            if ppsign==1,
                ensenya(['Loading data.']);
                [tagr,Ch,a,b,V,V2,nnr,RF]=CarregaVolums(folder,RF,vid,DT);
                [pos,posr,Q,Iryr,rr,rrr,sign,signals,scz,lim,signalsM,limmad]=ObtenirSenyalsNorm(G,folder, RF,DX,nnr,rs,a,b,V,DT,rp);
                if v(ssf)==12, signals=signals*1.04; end
                
                ensenya(['Events detection.']);
                win=ceil(85/DT); % ms (150 ms 180424)
                if v(ssf)==8, umbwW=0.2; if ssc==8, umbwW=0.45; end,end % 0.45 % if v(ssf)==10, umbwW=0.8; end
                if v(ssf)==11, if ssc==4, umbwW=0.45; end,end 
                if v(ssf)==10&&ssc==4, umbwW=0.6; end
                [events,signalsM,ev,eventsW]=EventsWavDetection(folder,RF,tww,Tww,win,signalsM,DT,rr,umbw,umbwW,nnr,pos,Iryr,V2, posr,scz,Ch,mam);
                umbwW=umbwWo; close all;
            else
                ensenya(['Loading events detection.']);
                RF=[RF '_Events']; tagr='ch01'; Ch=1;
                load([folder,'\',RF,'\volums.mat']); load([folder,'\',RF,'\signals.mat']);load([folder,'\',RF,'\ev_det.mat']);
                Q=G.Q{1}; pos=[Q(:,3),Q(:,2)]; rrr=length(pos); Iryr=G.Isum; rr=length(pos);     scz = get( 0, 'Screensize' ); nnr=size(V,3);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVENT PARAMETES %%%%%%%%%%%%%%%
            if ppdet==1,
                plt=0;
                [EV]=Imevents2MateventsTsort(events,pos,rr,Q); wts=round(100/DT);
                [EVcorr]=EventParameters(EV,folder,RF,nnr,signals,plt,DT,DX,wts,pos,V2,events);
                
                % EVcorr: [ROI ROIryr R/nR x y t sig(wt) evind indevcorr et amp AMP ror t2p dhm tau r2 bline];
                % ind:      1     2     3  4 5 6   7      8        9     10  11  12  13  14 15  16  17   18
                % EV =  [NumROIsenyal  NumROIryr  RyR/nRyR  x  y  Frame]
                %             1             2        3      4  5    6
            else
                dspko=dspk; wto=wt;
                load([folder,'\',RF,'\detectionParam.mat']);
                dspk=dspko; wt=wto;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUST & FILT %%%%%%%%%%%%%%%%%
            if ppefcl==1,
                ensenya('Filtering sparks.');
                factampo=factamp; ampFo=ampF; tauFo=tauF; t2pFo=t2pF;
                 if v(ssf)==8, factamp=factamp*0.5; end
                if v(ssf)==9, factamp=factamp*0.5; end
                if v(ssf)==10, if ssc>=10; factamp=factamp*0.5; end, end
                 if v(ssf)==12, factamp=factamp*0.5;  end, 
                 if v(ssf)==15&&ssc==3, tauF=[10 120]; t2pF=[3 24]; end
                % if v(ssf)==13, factamp=factamp*0.5;  end, 
                [EVcorr,blS,blSstd]=FilterEvents(folder,RF,ampF,tauF,rorF,t2pF,dhmF,r2F,EVcorr,ev,signals,eventsW,zzRyR,Iryr,pos,DX,DT,rr,factamp,nnr,limmad,lim);
                factamp=factampo; ampF=ampFo;tauF=tauFo; t2pF=t2pFo;
                guardaCSVsactivations(EVcorr,folder,RF,ampF,tauF,rorF,t2pF,dhmF,factamp);
               
                
                [EVfilt,Vdet]=ClustandRmvW(folder,RF,S,ss,factamp,DX,DT,pos,G,wt,dspk,rmw,EVcorr,EV,nnr,rr,Iryr,events,eventsW,signals,ampF,tauF,rorF,t2pF,dhmF,r2F,V2,scz,blS,blSstd,Ch,limmad,lim);
            else
                try,load([folder,'\',RF,'\ClustEv_' num2str(factamp) '.mat']); end%         try,load([folder,'\',RF,'ROIsAD_' num2str(factamp) '.mat']); end
                try,load([folder '\' RF '\RyRsandEv_' num2str(factamp) '.mat']);end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA PARAMETERS %%%%%%%%%%%%%
            if ppVol==1,
                if ~isempty(EVfilt),
                    L=[]; L2=[]; signals2=signals; signalsw=eventsW;  % EVfilto{:,24}=[EVfilt{:,24}]; EVfilt=EVfilto;
                    
                    ensenya('Trajectories');
                    for ii=1:size(EVfilt,1),EVfilt{ii,26}=1;end
                    [trajectories]=EventTrajectories(EVfilt,folder,RF,Iryr,pos);
                    distcov={}; eryr=[];% win=ceil(250/DT); % ms
                    [distcov, eryr]=distancesCovariance(EVfilt,Iryr,DX);
                    % eryr=[inR >1RyR 1stact numRyRs dt mean(D) max(D)];
                    %        1    2     3      4      5    6      7
                    % distcov=[D cent CovM [Sx,Sy] eigV sqrt(lambda) ve*lamb Iev];
                    %          1  2     3     4      5       6          7     8
                    save([folder,'\',RF,'\trajectories.mat'],'EVfilt','signals','trajectories','EVfilt','eryr','distcov');
                    
                    uu=unique([EVfilt{:,10}]); RF2='distcov';
                    GG=gray(256); RR=GG; RR(:,2:3)=0; GG(:,[1,3])=0; GG(end-2,:)=[1 0 0]; GG(end-1,:)=[1 0 1]; GG(end,:)=[1 1 1];
                    delta=.5; wst=round(65/DT);
                    if (exist([folder '\' RF '\' RF2],'dir')~=0),rmdir([folder '\' RF '\' RF2],'s'); end
                    if (exist([folder '\' RF '\' RF2],'dir')==0), mkdir([folder '\' RF '\' RF2]);end
                    figure(1), clf;
                    for ii=1:length(uu)
                        if size(distcov,2)>2,
                            if ~isempty(distcov{ii,8}),
                                ind=find([EVfilt{:,10}]==uu(ii));
                                mx=min([nnr,max([[EVfilt{ind,6}],[EVfilt{ind,6}]+round(wt/DT)])]);
                                tt=min([[EVfilt{ind,6}] [EVfilt{ind,6}]+round(wt/DT)]):mx;
                                IV=sum(V2(:,:,tt),3);
                                ax1=subplot(211); imagesc(distcov{ii,8}); axis image off; colormap(ax1,GG); title(['Event ' num2str(uu(ii))]);
                                ax2=subplot(212); imagesc(IV); axis image off; colormap(ax2,RR);
                                saveWysiwyg(1,[folder,'\',RF,'\',RF2,'\distcov' num2str(uu(ii)) '.png']);
                            end
                        end
                    end
                    close(1);
                    
                    ensenya('RyRs act parameters.')
                    [a,b]=size(Iryr);
                    % Ratio activació dels RyRs (RyR act / RyR total):
                    Iact=sum(Vdet,3); Iacc=Iact; Iacc(Iacc>0)=1;
                    IQ=zeros(size(Iryr)); for ii=1:length(pos), IQ(pos(ii,1),pos(ii,2))=1; end;
                    ratioRyRsact=numel(find((Iacc+IQ*2)>=3))/rr;
                    
                    % Freq activation
                    % Ev / RyR:
                    RyRsF=zeros(1,rr);
                    for ii=1:rr
                        num=Iact(pos(ii,1),pos(ii,2));
                        RyRsF(ii)=num/(nnr*DT); % activations / ms
                    end
                    RyRsFtotal=mean(RyRsF);
                    RyRsFmin=min(RyRsF); RyRsFmax=max(RyRsF);
                    
                    % Ev / (area*t) - Spark frequency:
                    sf=zeros(1,nnr);
                    uu=unique([EVfilt{:,10}]); uu(uu==0)=[];
                    for ii=1:length(uu), ind=find([EVfilt{:,10}]==uu(ii)); sf(EVfilt{ind(1),6})=1; end
                    spkF=numel(find(sf==1))/(rr*nnr*DT); % spk/(ms*RyRs)
                    % Restar waves al temps
                    ddw=0;
                    if ~isempty(find(eventsW==1)),
                        dtw=zeros(length(rr),1);
                        for ii=1:rr, dtw(ii)=length(find(eventsW(ii,:)==1));end
                        %max(dtw)
                        ddw=max(dtw);nnr2=nnr-ddw;
                        spkF=numel(find(sf==1))/(rr*nnr2*DT);
                    end
                    for ii=1:size(EVfilt,1), EVfilt{ii,26}=1; end; % salts --> preparar
                    
                    if ~isempty(uu),
                        %  RyR act / Jumps:
                        numRyRsEV=zeros(1,length(uu)); coloc=numRyRsEV; colocR=numRyRsEV;
                        angRyRs=zeros(length(uu),3);
                        aRyRs={};
                        P1RyR=0; P2RyR=0; P3RyR=0;P4RyR=0; Pmes4RyRs=0;
                        rd=ceil(0.33/DX);
                        for ii=1:length(uu)
                            indr=find([EVfilt{:,10}]==uu(ii));
                            RyRs=[EVfilt{indr,1}];
                            tryrs=[EVfilt{indr,6}];
                            % Num de RyRs:
                            aa=[];
                            numRyRsEV(ii)=length(unique(RyRs));
                            if length(unique(RyRs))==1, P1RyR=P1RyR+1; aa=1; end;
                            if length(unique(RyRs))==2, P2RyR=P2RyR+1; aa=2; end;
                            if length(unique(RyRs))==3, P3RyR=P3RyR+1; aa=3; end;
                            if length(unique(RyRs))==4, P4RyR=P4RyR+1; aa=4; end;
                            if length(unique(RyRs))>4, Pmes4RyRs=Pmes4RyRs+1; aa=5; end;
                            % Index event num de RyRs:
                            for jj=1:length(indr), EVfilt{indr(jj),26}=aa; end
                            % Salts activacions
                            angRyRs(ii,3)=length(RyRs); dt=[];
                            if length(RyRs)>1,
                                ryri=trajectories{ii,7};
                                traj=trajectories{ii,6};
                                activ=trajectories{ii,3};
                                try,dt=trajectories{ii,4};dt=[0,dt];end
                                
                                if ~isempty(dt),
                                    aa=zeros(1,length(RyRs)-1); vv=aa;
                                    amp=aa; ror=aa; t2p=aa; fdhm=aa; tau=aa;
                                    ct=1;
                                    for jj=1:length(dt)
                                        if jj==1,
                                            if length(dt)>1,ind=find(traj(:,1)==1);else, ind=[1:length(ryri)]; end
                                            if length(ind)>1, % hi ha més d'un RyR
                                                if length(ind)<3, % només dos ryrs
                                                    P1=[pos(ryri(ind(1)),1) ,pos(ryri(ind(1)),2)];
                                                    P2=[pos(ryri(ind(2)),1) ,pos(ryri(ind(2)),2)];
                                                    [a]=angol([0,1],[P1-P2]);
                                                    aa(ct)=a;
                                                    amp(ct)=EVfilt{activ(ind(1)),11};
                                                    ror(ct)=EVfilt{activ(ind(1)),13};
                                                    t2p(ct)=EVfilt{activ(ind(1)),14};
                                                    fdhm(ct)=EVfilt{activ(ind(1)),15};
                                                    tau(ct)=EVfilt{activ(ind(1)),16};
                                                    ct=ct+1;
                                                else % més de 2 RyRs!
                                                    d=squareform(pdist([pos(ryri(ind),1),pos(ryri(ind),2)]));
                                                    d(d==0)=100000;
                                                    for kk=1:length(ind)-1 % Busca les parelles més properes
                                                        mm=min(min(d)); [c,r]=find(d==mm);
                                                        P1=[pos(ryri(ind(r(1))),1) ,pos(ryri(ind(r(1))),2)];
                                                        P2=[pos(ryri(ind(r(2))),1) ,pos(ryri(ind(r(2))),2)];
                                                        [a]=angol([0,1],[P1-P2]);
                                                        aa(ct)=a;
                                                        amp(ct)=EVfilt{activ(ind(r(1))),11};
                                                        ror(ct)=EVfilt{activ(ind(r(1))),13};
                                                        t2p(ct)=EVfilt{activ(ind(r(1))),14};
                                                        fdhm(ct)=EVfilt{activ(ind(r(1))),15};
                                                        tau(ct)=EVfilt{activ(ind(r(1))),16};
                                                        ct=ct+1;
                                                        d(r(1),c(1))=100000;
                                                        d(r(2),c(2))=100000;
                                                    end
                                                end
                                            end
                                        else
                                            % Trobar en cada fila quin activa a quin
                                            ind=find(traj(:,jj)>0);
                                            % if ~iempty(ind),
                                            for kk=1:length(ind)
                                                ryr1=traj(ind(kk),jj);
                                                ryr2=ryri(ind(kk));
                                                ind2=find(ryri==ryr1);
                                                if ryr1~=ryr2,
                                                    P1=[pos(ryr1,1),pos(ryr1,2)];
                                                    P2=[pos(ryr2,1),pos(ryr2,2)];
                                                    [a]=angol([0,1],[P1-P2]);
                                                    ve=(sqrt(((P1(1)-P2(1))^2)+((P1(2)-P2(2))^2))*DX)/(dt(jj)*DT); % Velocitat event
                                                    aa(ct)=a; vv(ct)=ve; % pix/fr
                                                    amp(ct)=EVfilt{activ(ind2),11};
                                                    ror(ct)=EVfilt{activ(ind2),13};
                                                    t2p(ct)=EVfilt{activ(ind2),14};
                                                    fdhm(ct)=EVfilt{activ(ind2),15};
                                                    tau(ct)=EVfilt{activ(ind2),16};
                                                    ct=ct+1;
                                                end
                                            end
                                            % end
                                        end
                                        
                                    end
                                    vv2=zeros(length(vv),2);
                                    aa2=zeros(1,length(aa));
                                    for jj=1:length(aa)
                                        vv2(jj,1)=vv(jj);
                                        if aa(jj)>180,
                                            if (aa(jj)<315)&&(aa(jj)>225),
                                                aa2(jj)=2;
                                                vv2(jj,2)=2;
                                            else
                                                aa2(jj)=1;
                                                vv2(jj,2)=1;
                                            end
                                        else
                                            if (aa(jj)<135)&&(aa(jj)>45)
                                                aa2(jj)=2;
                                                vv2(jj,2)=2;
                                            else
                                                aa2(jj)=1;
                                                vv2(jj,2)=1;
                                            end
                                        end
                                    end
                                    angRyRs(ii,1)=length(find(aa2==1))/length(aa); % Proporció mateixa Zline
                                    angRyRs(ii,2)=length(find(aa2==2))/length(aa); % Proporció diff zline
                                    angRyRs(ii,4)=(length(find(aa2==1))/length(aa))*(length(RyRs)-1); % n mateixa Zline
                                    angRyRs(ii,5)=(length(find(aa2==2))/length(aa))*(length(RyRs)-1); % n diff zline
                                    angRyRs(ii,6)=mean(vv);
                                    aRyRs{ii,1}=aa; % Numero de RyRs
                                    aRyRs{ii,2}=length(tryrs)*DT; % Duració activacions
                                    aRyRs{ii,3}=angRyRs(ii,4); % SZ
                                    aRyRs{ii,4}=angRyRs(ii,5); % DZ
                                    aRyRs{ii,5}=vv2; % velocitat (dx/dt) i direcció salt
                                    aRyRs{ii,6}=[amp',vv2(:,2)]; % amp
                                    aRyRs{ii,7}=[ror',vv2(:,2)]; % ror
                                    aRyRs{ii,8}=[t2p',vv2(:,2)]; % t2p
                                    aRyRs{ii,9}=[fdhm',vv2(:,2)]; % fdhm
                                    aRyRs{ii,10}=[tau',vv2(:,2)]; % tau
                                end
                            else
                                aRyRs{ii,1}=aa;
                                aRyRs{ii,2}=[];
                                aRyRs{ii,3}=[];
                                aRyRs{ii,4}=[];
                                %                                 aRyRs{ii,6}=mean([EVfilt{indr,11}]);
                                %                                 aRyRs{ii,7}=mean([EVfilt{indr,13}]);
                                %                                 aRyRs{ii,8}=mean([EVfilt{indr,14}]);
                                %                                 aRyRs{ii,9}=mean([EVfilt{indr,15}]);
                                %                                 aRyRs{ii,10}=mean([EVfilt{indr,16}]);
                            end
                        end
                        % Probabilitats d'act de RyRs:
                        P1RyRt=P1RyR/length(uu);
                        P2RyRt=P2RyR/length(uu);
                        P3RyRt=P3RyR/length(uu);
                        P4RyRt=P4RyR/length(uu);
                        Pmes4RyRst=Pmes4RyRs/length(uu);
                        sumP=P1RyRt+P2RyRt+P3RyRt+P4RyRt+Pmes4RyRst;
                        
                        P1RyRti=P1RyR;%/length(uu);
                        P2RyRti=P2RyR;%/length(uu);
                        P3RyRti=P3RyR;%/length(uu);
                        P4RyRti=P4RyR;%/length(uu);
                        Pmes4RyRsti=Pmes4RyRs;%/length(uu);
                        
                        % Salts entre Z o en la mateixa:
                        ind=find(angRyRs(:,3)>0); indd=[]; try,indd=find(angRyRs(:,4)>0); end
                        sameZ=sum(angRyRs(ind,3));
                        
                        if ~isempty(indd),diffZ=sum(angRyRs(indd,4));
                            SZ=sameZ/(sameZ+diffZ);
                            DZ=diffZ/(sameZ+diffZ);
                        else
                            SZ=1; DZ=0;
                            diffZ=0;
                        end
                        
                        Rd=[]; % Diffusió
                        
                        % eryr=[inR >1RyR 1stact numRyRs dt mean(D) max(D)];
                        %        1    2     3      4      5    6      7
                        % distcov=[D cent CovM [Sx,Sy] eigV sqrt(lambda) ve*lamb Iev];
                        %          1  2     3     4      5       6          7     8
                        % EVfilt 23 (inbl) 24 (et) 25 (wav)
                        
                        % Primer s'activa el RyR:
                        acR=length(find(eryr(1,:)==1))/length(eryr); % ev. que s'inicien en RyRs
                        acnR=1-acR;
                        % Calcular percentatge d'events amb més d'un RyR implicat
                        acRR=length(find(eryr(2,:)==1))/length(eryr); % ev. que s'inicien en RyRs
                    else
                        
                        acR=0; acnR=0; acRR=0;
                        angRyRs=[]; aRyRs=[]; SZ=0; DZ=0;
                        RyRsFmin=0; RyRsFmax=0; spkF=0;  RyRsFtotal=0; RyRsF=zeros(1,rr); ddw=0;
                        P1RyRt=0; P2RyRt=0;P3RyRt=0; P4RyRt=0; Pmes4RyRst=0; coloc=[]; colocR=[];
                    end
                    
                    save([folder,'\',RF,'\Det3DEvents_parameters.mat'],'events','signals',...
                        'signals2','signalsw','L','L2','trajectories','EVfilt','acR','acRR',...
                        'distcov','eryr','P1RyRt','P2RyRt','P3RyRt','P4RyRt','Pmes4RyRst','coloc','colocR',...
                        'spkF','RyRsFtotal','RyRsFmin','RyRsFmax','RyRsF','aRyRs','angRyRs','SZ','DZ','ddw');
                else
                    trajectories={};
                    distcov={}; eryr=[];
                    save([folder,'\',RF,'\trajectories.mat'],'EVfilt','signals','trajectories','EVfilt','eryr','distcov');
                    L=[]; L2=[]; signals2=signals; signalsw=eventsW;
                    acR=0; acnR=0; acRR=0;
                    angRyRs=[]; aRyRs=[]; SZ=0; DZ=0;
                    RyRsFmin=0; RyRsFmax=0; spkF=0; RyRsFtotal=0; RyRsF=zeros(1,rr); ddw=0;
                    P1RyRt=0; P2RyRt=0;P3RyRt=0; P4RyRt=0;Pmes4RyRst=0; coloc=[]; colocR=[];
                    save([folder,'\',RF,'\Det3DEvents_parameters.mat'],'events','signals',...
                        'signals2','signalsw','L','L2','trajectories','EVfilt','acR','acRR',...
                        'distcov','eryr','P1RyRt','P2RyRt','P3RyRt','P4RyRt','Pmes4RyRst','coloc','colocR',...
                        'spkF','RyRsFtotal','RyRsFmin','RyRsFmax','RyRsF','aRyRs','angRyRs','SZ','DZ','ddw');
                end
                try,
                    % Histogrames:
                    % Iact=sum(Vdet,3);
                    nRyR=10; nframes=round(150/DT); ndist=2.5*2.5;
                    c2=[26 195 181]/255;c1=[.6 .1 .8];
                    figure(1);clf;
                    subplot(221)
                    hc=histc(eryr(4,:),[1:1:nRyR]);%n=length(size(eryr,2)); hc=hc/n;
                    h=bar([1:1:nRyR],hc);set(h,'FaceColor',c2,'FaceAlpha',.5); title('nº RyRs'); ylabel('nº Events'); xlabel('nº RyRs');
                    subplot(222)
                    hc=histc(eryr(5,:),[0:1:nframes]);%n=length(size(eryr,2)); hc=hc/n;
                    h=bar([0:1:nframes]*DT,hc);set(h,'FaceColor',c2,'FaceAlpha',.5); title('Peak time distance'); ylabel('nº Events');xlabel('dt (ms)');
                    subplot(223)
                    hc=histc(eryr(7,:),[0:.1:ndist]);%n=length(size(eryr,2)); hc=hc/n;
                    h=bar([0:.1:ndist],hc);set(h,'FaceColor',c2,'FaceAlpha',.5); title('Max distance'); ylabel('nº Events');xlabel('max dist (/um)');
                    subplot(224)
                    act=Iact(:);act(act==0)=[];
                    hc=histc(act,[1:1:rr]);%n=length(size(eryr,2)); hc=hc/n;
                    h=bar([1:1:rr],hc);set(h,'FaceColor',c2,'FaceAlpha',.5); title('Number of actiavions for each RyR'); ylabel('nº Events'); xlabel('Act/RyR');
                    saveWysiwyg(1,[folder,'\',RF,'\Event_hist.png']); close all;
                end
            end
            
            % end
        end
    end
end