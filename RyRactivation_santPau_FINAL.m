close all;
clear all;
clc;

%%%%%%%%%%%%%%%%% FOLDER MOUSE EXPERIMENTS
ensenya('Select the experiment folder');
FOLc = 'D:\Dades Lab\RyR-GFP-New'; % MOUSE FOLDER
% Folders normal strure:
%   > MICE (MVM_1)  -->  this one!
%       > CELLS (MVM_1_cell01) --> when running the program you can select this one
%           > SERIES (MVM_1_cell01_serie1)


% ACTIVATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Activation Simple %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
usuari=1;           % Program requires user interaction (1)
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

% Event grouping
wt=50;              % Temporal window (ms) wt=20; 181004 wt = 50;
ccorr=0.7;          % Correlation
dspk=2.5;%1;           % um

% Filtering parametes
zzRyR=0.33;               % Accepted RyRs mask (no RyRs in the borders accepted) in um
ampF=.36;               % F/Fo .6 1.2 181120 .1
factamp=6;      % Fact multiply bl
tauF=[6 200];          % ms 5 100
rorF=[0.005  1.5];     % 1/ms
t2pF=[3 100];          % ms 3 100
dhmF=[12 250];            % ms (6 190307)
r2F=.05;                % .29
plt=1;                  % Graphs sparks parameters
rmw=1;                  % Remove waves

% Reprocess
ppryr=1;            % RyRs detection
vid=1;              % Ca channel videos
ppsign=1;           % Signals
ppdet=1;          % Event paramerets
ppefcl=1;             % Filtering events and clustering
ppVol=1;            % Volume parameters
ppVdet=1;           % Detection 3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% CODE main
% RyRactivation_stpau.m & mass_sparks.m
% num_dx_totals.m
% trajectories_all.m


RFo=RF;
% Folders normal strure:
%   > MICE (MVM_1) 
%       > CELLS (MVM_1_cell01) --> This one now!
%           > SERIES (MVM_1_cell01_serie1)

if usuari==1,
    FOL=uigetdir([FOLc],'Select the cell folder'); S=dir([FOL '\*MVM*']);
else
    FOL=[FOLc '\' Sc(ssc).name]; % S=dir([FOL '\*RyR-GFP*']);
end

ensenya('Running the analysis for the series of the experiment');
for ss=1:size(S,1) % Series
    ensenya(S(ss).name);
    folder=[FOL '\' S(ss).name];%folder=[S{ss}];
    RF=RFo;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RyRs %%%%%%%%%%%%%%%%%%%%%%%%%%
    if ppryr==1,
        ensenya(['Detecting RyRs.']); thIo=thI; fro=fr;
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
             
        ensenya(['Events detection.']);
        win=ceil(85/DT); mam=0;
        [events,signalsM,ev,eventsW]=EventsWavDetection(folder,RF,tww,Tww,win,signalsM,DT,rr,umbw,umbwW,nnr,pos,Iryr,V2, posr,scz,Ch,mam);
        close all;
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
        ensenya('Filtering activations.');
        factampo=factamp; ampFo=ampF; tauFo=tauF; t2pFo=t2pF;
        
        [EVcorr,blS,blSstd]=FilterEvents(folder,RF,ampF,tauF,rorF,t2pF,dhmF,r2F,EVcorr,ev,signals,eventsW,zzRyR,Iryr,pos,DX,DT,rr,factamp,nnr,limmad,lim);
        factamp=factampo; ampF=ampFo;tauF=tauFo; t2pF=t2pFo;
        guardaCSVsactivations(EVcorr,folder,RF,ampF,tauF,rorF,t2pF,dhmF,factamp);
        
        [EVfilt,Vdet]=ClustandRmvW(folder,RF,S,ss,factamp,DX,DT,pos,G,wt,dspk,rmw,EVcorr,EV,nnr,rr,Iryr,events,eventsW,signals,ampF,tauF,rorF,t2pF,dhmF,r2F,V2,scz,blS,blSstd,Ch,limmad,lim);
    else
        try,load([folder,'\',RF,'\ClustEv_' num2str(factamp) '.mat']); end%         try,load([folder,'\',RF,'ROIsAD_' num2str(factamp) '.mat']); end
        try,load([folder '\' RF '\RyRsandEv_' num2str(factamp) '.mat']);end
    end
    % Results: Events_filt_... .cvs --> Accepted activations
    
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPARKS
    ensenya('Sparks Analysis.');
    freqi=0; densi=0; masst=[]; massspkt=[]; nspk2=0; sparkst=[]; freq=[]; dens=[];
    if ~isempty(EVfilt),
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
        signalso=signals; ctcell=0;
        
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
%             if vsf(sss)==12&&(ssc==10||ssc==9||ssc==11),
%                 SIGm=SIGm(1:round(length(SIGm)/2)); M=M(1:round(length(M)/2));
%                 Fo=median(SIGm-M*2); lim(kk)=Fo; Fo2=Fo;
%             end
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
            bg=min(min(tr0));
            tr0=tr0-bg;tr0(tr0<0)=0;
            tr0=(tr0)/(max(max(max(tr0))));
            [yy,xx]=find(tr0==max(max(tr0)));
            if(size(tr0,1)~=size(tr0,2)),extr=' candidate is too close to image limits.     ';else extr='.     ';end
            try
                % R=BlobRadiusC(tr0,[yy(1) xx(1) ce(3)]);
                R=DX*Blob4Diameters(tr0,1,.5,0)/2;
                I=sum(volum(:,:,ot2:ft2),3);
                I=I/(max(max(I))); I=I*255; Io=I;
                I(a1,b1-round(R/DX):b1+round(R/DX))=256;
                I(a1-round(R/DX):a1+round(R/DX),b1)=256;
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
            
%             if vsf(sss)==5&&ssc==4&&pos2(ii,4)==1||vsf(sss)==11&&ssc==4, SIGm2 = medfilt1(sig2,ww*6,'truncate');
%                 act=SIGm2(inmn:inmx);
%                 sigs=sig2(inmn:inmx);
%                 ssi=sigs-act;
%                 vals=ssi(ssi>0);
%             end
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
            sparks(ii,15)=ctcell;
        end
        
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
        save([folder '\' RF '\Sparks_parameters.mat'],'sparks','sparkst','massspk');
        freqi=nspk/(nnr*DT);
        densi=nspk/(rr*nnr*DT);
    end
    freq=[freq;freqi];
    dens=[dens;densi];
    try, save([folder '\' RF '\figs.mat'],'sparks','massspk','pos2','V','signals'); end
    ctcell=ctcell+1;
    
    T = array2table(sparks(:,[1:7,9:13]),'VariableNames',{'time_peak','amp','AMP','RoR','t2p','FDHM','tau','base_line','nRyRs','FWHM','mass','drug'});
    writetable(T ,[folder '\' RF  '\sparks_parameters.xlsx']);
    
end



