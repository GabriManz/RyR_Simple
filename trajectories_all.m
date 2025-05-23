% Trajectories drugs
clear all;

% Parameters
usuari=0;
factampt=6;
% Event grouping
wt1=350;              % Temporal window (ms) wt=20; 181004 wt = 50;
ccorr=0.5;          % Correlation
dspk1=2.5;%1;           % um

RF='ResRyR_12_2';

if usuari==1,
    FOL2=uigetdir([],'Select the main experiments folder');
else
    FOL2='D:\Dades Lab\RyR-GFP_Activation_Sparks'; mam=0;
    % FOL2='D:\Dades Lab\RyR-GFP_MVM_LS'; mam=0;
end
Sf=dir([FOL2 '\*RyR-GFP*']);

RFo=RF;
v=[1:16];
factamp=factampt(1); MMt=[]; ctf=1; FOLtt={};
% spk in waves
SPKint=[]; SPKnoint=[]; SPKint2=[]; SPKnoint2=[];
for ssf=1:length(v) % Mouse
    if usuari==1,
        FOLc=uigetdir([FOL2],'Select the main experiments folder'); Sc=dir([FOLc '\*RyR-GFP*']);
    else
        FOLc=[FOL2 '\' Sf(v(ssf)).name]; Sc=dir([FOLc '\*RyR-GFP*']);
    end
    % FOLc=[FOL2 '\' Sf(v(ssf)).name]; Sc=dir([FOLc '\*RyR-GFP*']);
    SPKin=zeros(size(Sc,1),7);
    SPKnoin=zeros(size(Sc,1),7);
    SPKin2=zeros(size(Sc,1),7);
    SPKnoin2=zeros(size(Sc,1),7);
    for ssc=1:size(Sc,1) % Cell
        if usuari==1,
            FOL=uigetdir([FOLc],'Select the main experiments folder'); S=dir([FOL '\*RyR-GFP*']);
        else
            FOL=[FOLc '\' Sc(ssc).name]; S=dir([FOL '\*RyR-GFP*']);
        end
        % FOL=[FOLc '\' Sc(ssc).name]; S=dir([FOL '\*RyR-GFP*']);
        S([S.isdir]==0)=[];
        
        ind1=[];
        for ss=1:size(S,1),
            indls=strfind(S(ss).name,'LS');
            imgse=strfind(S(ss).name,'Image');
            if ssf>5,sestim=strfind(S(ss).name,'Serie1');else, sestim=[];end
            if any([~isempty(indls),~isempty(imgse),~isempty(sestim)])
                ind1=[ind1;ss];
            end
        end
        if ~isempty(ind1), S(ind1)=[];end
        EVfilt3=[]; IRyRs=zeros(40,256); trajectoriest=[];
        
        % Carregar RyRs
        post=cell(1,size(S,1)); Iryrt=[];
        for ss=1:size(S,1)
            FOL=[FOLc '\' Sc(ssc).name]; % S=dir([FOL '\*RyR-GFP*']);
            folder=[FOL '\' S(ss).name];%folder=[S{ss}];
            RF=RFo;
            G=load([folder,'\',RF,'\finaldetection.mat']);
            Q=G.Q{1}; pos=[Q(:,3),Q(:,2)];
            DX=G.DX; DT=G.DT;
            post{1,ss}=pos;
            IQ=G.IQsum;
            if ~isempty(Iryrt), Iryrt=Iryrt+G.Inorm;else Iryrt=G.Inorm;end;
        end
        Iryrt=Iryrt/3;
        
        if size(S,1)>2,
            Ci=cell(1, size(S,1)-1);
            for ss=1:size(S,1)-1,Ci{ss}=union(post{1,ss},post{1,ss+1},'rows');end
            C=union(Ci{1},Ci{2},'rows'); if size(S,1)>3, C=union(C,Ci{3},'rows'); end
        else
            if size(S,1)==1,
                C=pos;
            else
                for ss=1:size(S,1)-1, C=union(post{1,ss},post{1,ss+1},'rows');end
            end
        end
        
        % Agrupa RyRs
        radi=ceil(0.33/DX); if DX>0.13, radi=ceil(0.5/DX); end
        D=pdist(C);
        Z =linkage(D);
        cc = cluster(Z,'cutoff',radi,'Criterion','distance');
        
        C=[C,cc];
        
        posR=zeros(max(cc),11);
        for ii=1:max(cc)
            ind=find(C(:,3)==ii);
            posR(ii,1)=round(mean(C(ind,1)));
            posR(ii,2)=round(mean(C(ind,2)));
        end
        [~,ia]=sort(posR(:,2));
        posR=posR(ia,:);
        
        %         I = zeros(40,256);
        %         for ii=1:size(posR,1), I(posR(ii,1),posR(ii,2))=1; end;
        %         imagesc(I*(max(max(G.Inorm)))+G.Inorm); axis image;
        radi=ceil(0.4/DX); if DX>0.13, radi=ceil(0.5/DX); end
        se=strel('disk',radi); se=double(se.Neighborhood);
        for ii=1:size(S,1),
            R=post{1,ii};
            for jj=1:size(R,1)
                I=zeros(size(G.Inorm));
                pp=R(jj,:);
                
                D=zeros(1,size(posR,1));
                for kk=1:size(posR,1),xx=posR(kk,1); yy=posR(kk,2); D(kk)=sqrt(((xx-pp(1))^2)+(yy-pp(2))^2); end
                
                [~,ia]= min(D);
                
                if posR(ia,ii+2)>0,
                    if posR(ia,ii+5)<1, posR(ia,ii+5)=jj;
                    else,  posR(ia,ii+8)=jj; end,
                else, posR(ia,ii+2)=jj; end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAVES!
        FOLtt{ctf,1}=Sc(ssc).name; ctf=ctf+1;
        % trajectories
        Att=[];Att2=[]; PLt=[]; dtwt=[];
        % Sparks iniciadors wave
        spkinwavet=[];
        spknoinwavet=[];
        spkinwavet2=[];
        spknoinwavet2=[];
        for ss=1:size(S,1) % Serie
            if usuari==1,
                FOL=uigetdir([FOLc],'Select the main experiments folder'); S=dir([FOL '\*RyR-GFP*']);
            else
                FOL=[FOLc '\' Sc(ssc).name]; % S=dir([FOL '\*RyR-GFP*']);
            end
            
            folder=[FOL '\' S(ss).name];%folder=[S{ss}];
            RF=RFo;
            G=load([folder,'\',RF,'\finaldetection.mat']);
            DX=G.DX; DT=G.DT;
            RF=[RF '_Events']; tagr='ch01'; Ch=1;
            load([folder,'\',RF,'\volums.mat']); load([folder,'\',RF,'\signals.mat']);load([folder,'\',RF,'\ev_det.mat']);
            Q=G.Q{1}; pos=[Q(:,3),Q(:,2)]; rrr=length(pos); Iryr=G.Isum; rr=length(pos);     scz = get( 0, 'Screensize' ); nnr=size(V,3);
            load([folder,'\',RF,'\detectionParam.mat']);
            try,load([folder,'\',RF,'\ClustEv_' num2str(factamp) '.mat']); end
            try,load([folder '\' RF '\RyRsandEv_' num2str(factamp) '.mat']);end
            EVfilto=EVfilt; signalso=signals;
            try, load([folder '\' RF '\figs.mat']); end % sparks
            signals=signalso;
            
            % drug:
            ind=find(folder=='_'); ind1=ind(end-2);ind2=ind(end-1); drug=folder(ind1+1:ind2-1); indx=0;
            if strfind(drug,'RyR')==1, indx=1;  end % CON
            if strfind(drug,'CIL')==1, indx=2;  end, if strfind(drug,'ISO')==1, indx=3; end %% CIL CIL+CPA
            if strfind(drug,'RO')==1, indx=4; end, if strfind(drug,'FENO')==1, indx=5; end % RO RO+CPA
            if strfind(drug,'RO+CIL')==1, indx=6; end, if strfind(drug,'CGS')==1, indx=7; end % RO+CIL CGS
            if indx==0, ind1=ind(end-3);ind2=ind(end-2); drug=folder(ind1+1:ind2-1); % hi ha algun comentari
                if ~isempty(strfind(drug,'RyR')==1)||~isempty(strfind(drug,'pl2')==1), indx=1; end % CON
                if strfind(drug,'CIL')==1, indx=2; end,  if strfind(drug,'CIL+CPA')==1, indx=3; end,  if strfind(drug,'RO')==1, indx=4; end, if strfind(drug,'RO+CPA')==1, indx=5; end, if strfind(drug,'RO+CIL')==1, indx=6; end, if strfind(drug,'CGS')==1, indx=7; end % CGS
            end
            
            dtwave=[];
            if any(any(eventsW>0)),
                % Efilt --> ajuntar info sparks i waves
                %       --> ajuntar events diff cells
                
                eventsW=conv2(eventsW,ones(1,round(50/DT)),'same'); eventsW(eventsW>0)=1;
                EVwave=cell(1,size(EVfilt,2));ct=1;
                for ii=1:rr
                    sig=eventsW(ii,:);
                    if any(sig>0),
                        bw=bwlabel(sig); % events separats
                        for jj=1:max(bw)
                            EVwave{ct,1}=ii; % ROI RyRs
                            EVwave{ct,2}=Q(ii,1); % tag RyRs
                            EVwave{ct,3}=1;
                            EVwave{ct,4}=pos(ii,1); % x
                            EVwave{ct,5}=pos(ii,2); % y
                            EVwave{ct,6}=find(bw==jj,1,'first'); % t
                            EVwave{ct,7}=find(bw==jj,1,'last');
                            % EVwave{ct,7}=sig();
                            ct=ct+1;
                        end
                    end
                end
                EV=[EVwave{:,1:7}]; EV=reshape(EV,[],7);
                
                % Calcula paràmetres wave
                %%% signals2=eventsW.*signals; % imagesc(signals2);
                plt=0;
                [EVw]=Imevents2MateventsTsort(eventsW,pos,rr,Q); wts=round(100/DT);
                [EVcorrwaves]=EventParameters(EVw,folder,RF,nnr,signals,plt,DT,DX,wts,pos,V2,events);
                
               
                
                % clusteritzar waves
                [EVwave,et,Iryr]=eventsRyRsClust4(wt1,DT,DX,dspk1,signals,Iryr,rr,nnr,EV,EVwave);
                mm=0; if ~isempty(EVfilt),mm=max([EVfilt{:,10}]); end, if isempty(mm),mm=0; end
                for ii=1:size(EVwave,1), EVwave{ii,10}=EVwave{ii,10}+mm; end
                EVfilt2=[EVfilt;EVwave];
                
                % Sparks before waves
                inds=[];
                etw=unique([EVwave{:,10}]);
                for ii=etw
                    ind=find([EVwave{:,10}]==ii);
                    ttw=min([EVwave{ind,6}]); % temps iniciador wave
                    % Buscar el nombre més gran dels tt<0
                    ttspk=[sparks(:,1)]-ttw; 
                    ttspk(ttspk>0)=[]; % es queda amb els t<0 (abans de la wave)
                    [~,c]=min(abs(ttspk)); % el més proper a la wave
                    % Guardar sparks (abans de la waves 150 ms)
                    nryrs=sparks(c,10);
                    if abs(ttspk(c))<round(300/DT),
                        inds=[inds;[c,nryrs]];
                    end
                end
                if ~isempty(inds),
                    vv=1:size(sparks,1); 
                    for jj=1:length(inds(:,1)),vv(vv==inds(jj,1))=[];end
                    if ~isempty(inds),
                        spkinwave=mean(sparks(inds(:,1),[2,4:7,10]),1); % iniciadors wave
                        spknoinwave=mean(sparks(vv,[2,4:7,10]),1); % no iniciadors wave
                    end
                    % Mirar si nombre de RyRs és diferent
                    % Mirar si amp o FDHM és més gran
                    % per cell:
                    spkinwavet=[spkinwavet;spkinwave];
                    spknoinwavet=[spknoinwavet;spknoinwave];
                end
                
                % Upraise calcium:
                [EVfilt]=UpRiseCalcium(EVfilt2, DT,signals);
                
                % dt inici ryr1 - ryrn wave (ms)
                mm=unique([EVwave{:,10}]);
                dtwave=zeros(1,length(mm));
                for ii=1:length(mm)
                    ind=find([EVwave{:,10}]==mm(ii));
                    if length(ind)>rr*0.15,%*0.5,
                        tt=[EVwave{ind,6}];
                        dtwave(ii)=[max(tt)-min(tt)]*DT/rr; % ms/nryrs
                    end
                end
                dtwave(dtwave==0)=[]; dtwave=mean(dtwave);
                
                IQ=zeros(size(Iryr));
                for ii=1:size(Q),IQ(Q(ii,3),Q(ii,2))=1;end
                
                delta=.5; wst=round(100/DT);
                % CLUSTERS AND ACTIVATIONS
                Qr=G.Q; Qr=Qr{1};
                ind=unique([EVfilt{:,1}]);
                labb=unique([EVfilt{:,10}]);
                cmclust=superjet(length(labb),'lines');
                if labb==1, cmclust=[0.9290, 0.6940, 0.1250];end
                tts=1:nnr;
                figure(1),clf;set(gcf,'Position',scz);
                for jj=1:rr
                    s=signals(jj,:)+jj*delta;
                    plot(tts*DT,s,'k'); hold on;
                    inde=find([EVfilt{:,1}]==jj);
                    if ~isempty(inde),
                        for kk=1:length(inde)
                            c=EVfilt{inde(kk),6};
                            lab=EVfilt{inde(kk),10};
                            indl=find(labb==lab);
                            tt=[max([1 c-5]):min([nnr c+wst+5])];
                            % Cluster
                            plot(tt*DT,s(tt),'Color',cmclust(indl,:),'LineWidth',1.6);
                            
                            % Act / Diff
                            cool='.r'; %if [EVfilt{inde(kk),25}]==1, cool='.y'; end % ETIQUETA wave
                            scatter(c*DT,s(c)+.2,cool);
                            scatter(EVfilt{inde(kk),23}*DT,s(EVfilt{inde(kk),23})+.2,15,'k','filled');
                            %pause;
                            mm=mean(s);
                            text((tts(1)*DT)-20,mm,num2str(Qr(EVfilt{inde(kk),1},1)),'Color','k');
                        end
                    end
                end
                xlabel('t (ms)'); set(gca,'ytick',[]);
                saveWysiwyg(1,[folder,'\',RF,'\Signals_clust_AD_Waves_' S(ss).name '_' num2str(factamp) '.png']);
                
            end
            
            % Trajectories
            if ~isempty(EVfilt),
                ensenya('Trajectories');
                for ii=1:size(EVfilt,1),EVfilt{ii,26}=1;end
                [trajectories]=EventTrajectories(EVfilt,folder,RF,Iryr,pos);
                
                
                % Adjency matrix per a cada event
                ml=1;
                rr2=size(posR,1); PL=zeros(1,rr2);
                At=zeros(rr2); At2=At;
                if size(trajectories,2)<4,
                    for ii=1:size(trajectories,1)
                        A=zeros(rr2); % activacions, central és el primer activat
                        A2=zeros(rr2); % activats sols
                        uu=unique([EVfilt{:,10}]);et=uu(ii);
                        ind1=find([EVfilt{:,10}]==et);
                        ind=find(posR(:,ss+2)==EVfilt{ind1,1});
                        A2(ind,ind)=1; % Activador inicial
                        trajectories{ii,8}=A; % Mutiple clust act
                        trajectories{ii,9}=A2; % single clust act
                        At=At+A2; % all events
                        At2=At2+A2;
                    end
                else
                    for ii=1:size(trajectories,1)
                        A=zeros(rr2); % activacions, central és el primer activat
                        A2=zeros(rr2); % activats sols
                        A3=zeros(rr2);
                        M=trajectories{ii,6}; % Trajectoria activacions RyRs - intervals temps
                        vv1=trajectories{ii,7}; % RyRs en ordre
                        vv=[];
                        for jj=1:size(vv1,2), ind=find(posR(:,ss+2)==vv1(jj));
                            if isempty(ind),ind=find(posR(:,ss+5)==vv1(jj)); end
                            if isempty(ind),ind=find(posR(:,ss+8)==vv1(jj)); end
                            vv=[vv;ind];
                        end
                        if isempty(M), % només s'ha activat un/dos
                            if isempty(vv), % un RyRs
                                uu=unique([EVfilt{:,10}]);
                                et=uu(ii);
                                ind1=find([EVfilt{:,10}]==et);
                                ind=find(posR(:,ss+2)==EVfilt{ind1,1});
                                A(ind,ind)=1; % Activador inicial
                               % A3(ind,ind)=1;
                            else % més d'un RyR co-activat
                                for kk=1:length(vv)
                                    A(vv(kk),vv(kk))=1; % Activador inicial
                                    A3(vv(kk),vv(kk))=1;
                                    if PL(vv(kk))<length(vv), PL(vv(kk))=length(vv); end
                                end
                            end
                        else % Activació de diferents clusters
                            ind=find(M(:,1)==1);
                            [~,cc]=max(PL(ind));
                            for kk=1:length(ind),
                                A(vv(ind(kk)),vv(ind(kk)))=1; % Activador inicial
                                A3(vv(ind(kk)),vv(ind(kk)))=1;
                                indrr=vv(ind(kk));
                                if kk==cc,
                                    if PL(vv(ind(kk)))<size(M,2), PL(vv(ind(kk)))=size(M,2); end % Path length
                                end
                            end
                            for kk=2:size(M,2)
                                ind=find(M(:,kk)>0); ind2=zeros(size(ind,1),1);
                                for jj=1:length(ind),ind2(jj)=find(vv1==M(ind(jj),kk),1,'first'); end
                                A(vv(ind2),vv(ind))=1; % Activacions de cada RyR (ex. 1 act 2, 2 act 4)
                                A3(indrr,vv(ind))=1; % Activacions del RyR activador (ex. 1 act 2, 1 act 4 (a través de 3, però el posa en el primer activador)
                            end
                            % max length
                            ml=[ml;size(M,2)];
                        end
                        trajectories{ii,8}=A; % Mutiple clust act
                        trajectories{ii,9}=A2; % single clust act
                        At=At+A+A2; % all events
                        At2=At2+A3+A2;
                    end
                end
                % imagesc(At2); axis square;
                
                %max length
                ml=max(PL); % máx nombre de t
                
                pos2=[pos,ones(length(pos),1)*ss];
            end
            
            if ~isempty(EVfilt)
                % canviar et
                if ~isempty(EVfilt3),
                    mm=max([EVfilt3{:,10}]);
                    for ii=1:size(EVfilt,1), EVfilt{ii,10}=EVfilt{ii,10}+mm; end
                    if size(EVfilt,2)<size(EVfilt3,2), EVfilt{size(EVfilt3,1),size(EVfilt3,2)}=0;end
                    EVfilt3=[EVfilt3;EVfilt]; % acumula dades
                    trajectoriest=[trajectoriest;trajectories];
                    if isempty(Att), Att=At; else, Att=Att+At;end
                    if isempty(Att2), Att2=At2; else; Att2=Att2+At2;end
                    if isempty(PLt), PLt=PL; else PLt=PLt+PL; end
                    if isempty(dtwt), dtwt=dtwave; else, dtwt=[dtwt; dtwave]; end
                else
                    EVfilt3=EVfilt; 
                    trajectoriest=[trajectories];
                    if isempty(Att), Att=At; else, Att=Att+At;end
                    if isempty(Att2), Att2=At2; else; Att2=Att2+At2;end
                    if isempty(PLt), PLt=PL; else PLt=PLt+PL; end
                    if isempty(dtwt), dtwt=dtwave; else, dtwt=[dtwt; dtwave]; end
                end
            end
        end
        if ~isempty(spkinwavet),
            SPKin(ssc,:)=[mean(spkinwavet,1),indx];
            SPKnoin(ssc,:)=[mean(spknoinwavet,1),indx];
        else
            SPKin(ssc,7)=indx;
            SPKnoin(ssc,7)=indx;
        end
        
        
        
        if ~isempty(Att),
        figure(1),clf;
        subplot(121); imagesc(Att);axis square; title('Adjacency Matrix RyRs cl. to RyRs cl.');
        subplot(122); imagesc(Att2);axis square; title('Adjacency Matrix 1st RyRs cl. activated');
        saveWysiwyg(1,[FOL '\trajectories_AdjacencyMatrix.png']);
        
        % max path length (de les series)
        mpl=max(PLt)/rr2;
        
        % RyRs iniciadors
        for ii=1:size(Att2,1)
            Ra(ii)=Att2(ii,ii);
        end
        Ract=numel(find(Ra>0));
        Ractt=Ract/rr2; % percentatge de RyRs activadors
        Ra(Ra==0)=[];
        nactRa=mean(Ra); % Nombre activacions RyRs activadors
        
        % RyRs activats
        ryrs=sum(Att,1); % percentatge de RyRs activats en totes les seqüencies
        ryrst=numel(find(ryrs>0))/rr2;
        
        % Max node (max número d'activacions d'un RyR)
        mnode=max(max(Att));
        
        % Max node & max path length
        if ~isempty(PLt),mNP=max(min([max(Att)', PLt']'));else, mNP=0; end
        
        % path length / max node
        % if ~isempty(PLt),vv= PLt'./max(Att)'; vv(isinf(vv))=[];rNP=max(vv);rrNP=mean(vv(vv>0));else, rNP=0;rrNP=0; end
        if ~isempty(PLt),vv= PLt'.*max(Att)'; vv(isinf(vv))=[];rNP=max(vv);
            if any(vv)>0, rrNP=mean(vv(vv>0)); else, rrNP=0;end
        else, rNP=0;rrNP=0; 
        end
        
        
        % dt waves
        if ~isempty(dtwt) dtw=mean(dtwt); else, dtw=0; end
        
        cmap=gray(256); cmap(:,[1,3])=0;
        cmap2=jet(max(ryrs)+1);
        cmap3=[1 1 1];
        
        figure(1); clf;
        ax1=subplot(211); imagesc(Iryrt);hold on; colormap(ax1,cmap); axis image off;
        subplot(212);
        imagesc(ones(size(Iryrt))); colormap(cmap3); 
        hold on;
        for ii=1:size(posR,1)
            ind=find(Att(ii,:)>0);
            if ~isempty(ind),
            for jj=1:length(ind)
                plot([posR(ind(jj),2),posR(ii,2)],[posR(ind(jj),1),posR(ii,1)],'Color',[0.5 0.5 0.5],'LineWidth',1.5)
            end
            end
        end
         for ii=1:size(posR,1)
            %if any(Att(ii,:)>0)
            scatter(posR(ii,2),posR(ii,1),80,'MarkerEdgeColor',cmap2(ryrs(ii)+1,:),...
              'MarkerFaceColor',[1 1 1],'LineWidth',2);
            % end
        end
        set(gcf,'position',[1 41 1920/2.5 964/3]);axis image off;
        saveWysiwyg(gcf,([FOL '\connectivity.png']));
        
        else
            ryrst=0; Ractt=0; mpl=0; mnode=0; nactRa=0; mNP=0;dtw=0; rNP=0;rrNP=0;
        end
        
         
        
        % Traj summary
        MM=[ryrst,Ractt,nactRa,mpl,mnode,mNP,dtw,rr2,indx,rNP,rrNP];
        
        save([FOL '\trajectories.mat'],'EVfilt3','Att','Att2','posR','trajectoriest','MM','indx');
        MMt=[MMt;MM];
    end
    SPKint=[SPKint;SPKin];
    SPKnoint=[SPKnoint;SPKnoin];
end

% spk in wave
icon=find(SPKint(:,7)==1);
ifeno=find(SPKint(:,7)==5);
iiso=find(SPKint(:,7)==3);
Tcon = array2table(SPKint(icon,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
Tfeno = array2table(SPKint(ifeno,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
Tiso = array2table(SPKint(iiso,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
writetable(Tcon,[FOL2 '\spk_inicio_waves_CON.xlsx']);
writetable(Tfeno,[FOL2 '\spk_inicio_waves_FENO.xlsx']);
writetable(Tiso,[FOL2 '\spk_inicio_waves_ISO.xlsx']);
icon=find(SPKnoint(:,7)==1);
ifeno=find(SPKnoint(:,7)==5);
iiso=find(SPKnoint(:,7)==3);
Tcon = array2table(SPKnoint(icon,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
Tfeno = array2table(SPKnoint(ifeno,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
Tiso = array2table(SPKnoint(iiso,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
writetable(Tcon,[FOL2 '\spk_no_inicio_waves_CON.xlsx']);
writetable(Tfeno,[FOL2 '\spk_no_inicio_waves_FENO.xlsx']);
writetable(Tiso,[FOL2 '\spk_no_inicio_waves_ISO.xlsx']);


% recuperar totes les dades FOL
MMt=[];
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
        load([FOL '\trajectories.mat']);
        MMt=[MMt;MM];
    end
end

% MMt = [ryrst,Ractt,nactRa,mpl,mnode,mNP,dtw,RyRs,indx,max(PL/nNodes),<PL/nNodes>]

TT=table(FOLtt,MMt);
writetable(TT,[FOL2 '\Trajectories.xlsx']);

icon=find(MMt(:,9)==1);
iro=find(MMt(:,9)==4);
icil=find(MMt(:,9)==2);
iiso=find(MMt(:,9)==3);
ifeno=find(MMt(:,9)==5);

% scatter(MMt(:,4),MMt(:,3))


c33=[255, 204, 0]/255;
c55=[0, 134, 179]/255;
c55=[205, 237, 76]/255;
c11=[115, 0, 153]/255;
c44=[153, 0, 0]/255;
c22=[102, 153, 255]/255;
FS=14;

figure(1);clf;
v=[1,2,3,4,5,11];
noms={'RyRs activated','RyRs initiators','Mean activations RyRs initiators',...
    'Max path length','max #activations per RyRs','PathLength·#RyRactivations'};
uni={'#RyRsactivated/#RyRstotal','#RyRsinitiators/#RyRstotal','<activations>',...
    '#RyRsEvent/#RyRstotal','n','#RyRsPath·#activarions'};
for ii=1:length(v)
    x1=MMt(icon,v(ii));
    x2=MMt(icil,v(ii));
    x4=MMt(iro,v(ii));
    x5=MMt(ifeno,v(ii)); 
    x3=MMt(iiso,v(ii));
    subplot(2,3,ii)
    CompareDistributions({x1,x4,x2,x5,x3},'tst','ranksum','lab',{'CON','RO','CIL','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{ii},'sig','*', 'uni',uni{ii},'fig',0,...
    'fs',FS,'dc',[c11;c44;c22;c33;c55],'wi',0.5); axis square;
end
set(gcf,'position',[1 41 1920/2 964/1.5]);
saveWysiwyg(gcf,([FOL2 '\trajectories.png']));

figure(1);clf;
v=[1,2,3,4,5,11];
noms={'RyRs activated','RyRs initiators','Mean activations RyRs initiators',...
    'Max path length','max #activations per RyRs','PathLength·#RyRactivations'};
uni={'#RyRsactivated/#RyRstotal','#RyRsinitiators/#RyRstotal','<activations>',...
    '#RyRsEvent/#RyRstotal','n','#RyRsPath·#activarions'};
for ii=1:length(v)
    x1=MMt(icon,v(ii));
    x2=MMt(icil,v(ii));
    x4=MMt(iro,v(ii));
    x5=MMt(ifeno,v(ii)); 
    x3=MMt(iiso,v(ii));
    subplot(2,3,ii)
    CompareDistributions({x1,x4,x2,x3},'tst','ranksum','lab',{'CON','RO','CIL','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{ii},'sig','*', 'uni',uni{ii},'fig',0,...
    'fs',FS,'dc',[c11;c44;c22;c55],'wi',0.5); axis square;
end
set(gcf,'position',[1 41 1920/2 964/1.5]);
saveWysiwyg(gcf,([FOL2 '\trajectories_ISO.png']));



% waves sense CON
figure(1);clf;
noms={'fract. w / waves'}; uni={'n waves/n total'};
subplot(2,3,1)
x2=find(MMt(icil,7)>0); x2=length(x2)/length(icil);
x4=find(MMt(iro,7)>0); x4=length(x4)/length(iro);
x5=find(MMt(ifeno,7)>0); x5=length(x5)/length(ifeno);
x3=find(MMt(iiso,7)>0); x3=length(x3)/length(iiso);
CompareDistributions({x4,x2,x5,x3},'tst','ranksum','lab',{'RO','CIL','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{1},'sig','*', 'uni',uni{1},'fig',0,...
    'fs',FS,'dc',[c44;c22;c33;c55],'wi',0.5); axis square;
noms={'clusters in waves / time'}; uni={'maxPL/\DeltatWaves'};
% x1=(MMt(icon,4).*MMt(icon,8))./MMt(icon,7);x1(isnan(x1))=[]; x1(isinf(x1))=[]; 
x2=(MMt(icil,4).*MMt(icil,8))./MMt(icil,7); x2(isnan(x2))=[]; x2(isinf(x2))=[];
x4=(MMt(iro,4).*MMt(iro,8))./MMt(iro,7); x4(isnan(x4))=[]; x4(isinf(x4))=[];
x5=(MMt(ifeno,4).*MMt(ifeno,8))./MMt(ifeno,7); x5(isnan(x5))=[]; x5(isinf(x5))=[];
x3=(MMt(iiso,4).*MMt(iiso,8))./MMt(iiso,7); x3(isnan(x3))=[]; x3(isinf(x3))=[];
subplot(2,3,2)
CompareDistributions({x4,x2,x5,x3},'tst','ranksum','lab',{'RO','CIL','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{ii},'sig','*', 'uni',uni{ii},'fig',0,...
    'fs',FS,'dc',[c44;c22;c33;c55],'wi',0.5); axis square;
v=[7]; noms={'\Deltat waves / #RyRs'}; uni={'ms/clusters'};
    % x2=MMt(icil,v(ii)).*MMt(icil,v(ii)+1); x4=MMt(iro,v(ii)).*MMt(iro,v(ii)+1);  x5=MMt(ifeno,v(ii)).*MMt(ifeno,v(ii)+1);  x3=MMt(iiso,v(ii)).*MMt(iiso,v(ii)+1);
    x2=MMt(icil,v(ii));
    x4=MMt(iro,v(ii));
    x5=MMt(ifeno,v(ii)); 
    x3=MMt(iiso,v(ii));
    x1(x1==0)=[];x2(x2==0)=[];x3(x3==0)=[];x4(x4==0)=[];x5(x5==0)=[];
    subplot(2,3,3)
    CompareDistributions({x4,x2,x5,x3},'tst','ranksum','lab',{'RO','CIL','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{ii},'sig','*', 'uni',uni{ii},'fig',0,...
    'fs',FS,'dc',[c44;c22;c33;c55],'wi',0.5); axis square;
set(gcf,'position',[1 41 1920/2 964/1.5]);
saveWysiwyg(gcf,([FOL2 '\trajectories_waves.png']));

% waves sense CON
figure(1);clf;
noms={'fract. w / waves'}; uni={'n waves/n total'};
subplot(2,3,1)
x2=find(MMt(icil,7)>0); x2=length(x2)/length(icil);
x4=find(MMt(iro,7)>0); x4=length(x4)/length(iro);
x5=find(MMt(ifeno,7)>0); x5=length(x5)/length(ifeno);
x3=find(MMt(iiso,7)>0); x3=length(x3)/length(iiso);
CompareDistributions({x4,x2,x3},'tst','ranksum','lab',{'RO','CIL','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{1},'sig','*', 'uni',uni{1},'fig',0,...
    'fs',FS,'dc',[c44;c22;c55],'wi',0.5); axis square;
noms={'clusters in waves / time'}; uni={'maxPL/\DeltatWaves'};
% x1=(MMt(icon,4).*MMt(icon,8))./MMt(icon,7);x1(isnan(x1))=[]; x1(isinf(x1))=[]; 
x2=(MMt(icil,4).*MMt(icil,8))./MMt(icil,7); x2(isnan(x2))=[]; x2(isinf(x2))=[];
x4=(MMt(iro,4).*MMt(iro,8))./MMt(iro,7); x4(isnan(x4))=[]; x4(isinf(x4))=[];
x5=(MMt(ifeno,4).*MMt(ifeno,8))./MMt(ifeno,7); x5(isnan(x5))=[]; x5(isinf(x5))=[];
x3=(MMt(iiso,4).*MMt(iiso,8))./MMt(iiso,7); x3(isnan(x3))=[]; x3(isinf(x3))=[];
subplot(2,3,2)
CompareDistributions({x4,x2,x3},'tst','ranksum','lab',{'RO','CIL','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{1},'sig','*', 'uni',uni{1},'fig',0,...
    'fs',FS,'dc',[c44;c22;c55],'wi',0.5); axis square;
v=[7]; noms={'\Deltat waves / #RyRs'}; uni={'ms/clusters'};
    % x2=MMt(icil,v(ii)).*MMt(icil,v(ii)+1); x4=MMt(iro,v(ii)).*MMt(iro,v(ii)+1);  x5=MMt(ifeno,v(ii)).*MMt(ifeno,v(ii)+1);  x3=MMt(iiso,v(ii)).*MMt(iiso,v(ii)+1);
    x2=MMt(icil,v(ii));
    x4=MMt(iro,v(ii));
    x5=MMt(ifeno,v(ii)); 
    x3=MMt(iiso,v(ii));
    x1(x1==0)=[];x2(x2==0)=[];x3(x3==0)=[];x4(x4==0)=[];x5(x5==0)=[];
    subplot(2,3,3)
    CompareDistributions({x4,x2,x3},'tst','ranksum','lab',{'RO','CIL','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{1},'sig','*', 'uni',uni{1},'fig',0,...
    'fs',FS,'dc',[c44;c22;c55],'wi',0.5); axis square;
set(gcf,'position',[1 41 1920/2 964/1.5]);
saveWysiwyg(gcf,([FOL2 '\trajectories_waves_ISO.png']));

mean(x4)
std(x4)/sqrt(length(x4))
mean(x2)
std(x2)/sqrt(length(x2))
mean(x3)
std(x3)/sqrt(length(x3))
mean(x5)
std(x5)/sqrt(length(x5))

figure(1);clf;
v=[7];
noms={'\Deltat waves / #RyRs'};
uni={'ms/clusters'};
for ii=1:length(v)
    x1=MMt(icon,v(ii));
    x2=MMt(icil,v(ii));
    x4=MMt(iro,v(ii));
    x5=MMt(ifeno,v(ii)); 
    x3=MMt(iiso,v(ii));
    x1(x1==0)=[];x2(x2==0)=[];x3(x3==0)=[];x4(x4==0)=[];x5(x5==0)=[];
    % x1(1)=[];
    subplot(2,3,ii)
    CompareDistributions({x1,x4,x2,x3},'tst','ranksum','lab',{'CON','RO','CIL','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{ii},'sig','*', 'uni',uni{ii},'fig',0,...
    'fs',FS,'dc',[c11;c44;c22;c55],'wi',0.5); axis square;
end
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOL2 '\trajectories_dt_waves.png']));


% MMt = [ryrst,Ractt,nactRa,mpl,mnode,mNP,dtw,RyRs,indx,max(PL/nNodes),<PL/nNodes>]
% Fraccion de cèl·lulas que alcanzan el 10, 20, 50 y 90 % de activacion
figure(1);clf;
noms={'\Deltat waves / #RyRs'};
uni={'ms/clusters'};
bins=[0,.25,.5,0.75,1]*100; y=[];
indd={icon,iro,icil,iiso};
col=[c11;c44;c22;c55];hold on;
for ii=1:4
    ind=indd{ii};
    x1=MMt(ind,1);
   [h,c]=hist(x1*100,bins); h=h/sum(h);
   plot(c,h,'Color',col(ii,:),'LineWidth',2); 
   y=[y;h];
   
end
set(gcf,'position',[1 41 1920/2.5 964/2]);

c ={'10%','25%','50%','75%','100%'};
b=bar(y');
for ii=1:4
b(ii).EdgeColor=col(ii,:);
b(ii).FaceColor=col(ii,:);end
set(gca,'xticklabel',c)
xlabel('%RyRs activated per cell'); ylabel('fraction');
box off; % axis square;
set(gca,'FontSize',16);
saveWysiwyg(gcf,([FOL2 '\trajectories_RyRsact_ISO.png']));

% Max fracción de clusters activados / fracción de celulas con esta fracción.
fractcells=zeros(1,5);
for ii=1:4, ind=find(MMt(indd{ii},1)>0); fractcells(ii)=length(ind)/length(indd{ii}); end
figure(1);clf;
v=[1];
noms={'#Actiated RyRs / frac. active cells'};
uni={'#clust/fractcells'};
for ii=1:length(v)
    x1=MMt(icon,v(ii))*fractcells(1);
    x2=MMt(icil,v(ii))*fractcells(3);
    x4=MMt(iro,v(ii))*fractcells(2);
    x5=MMt(ifeno,v(ii))*fractcells(4); 
    x3=MMt(iiso,v(ii))*fractcells(5);
    
    %subplot(2,3,ii)
    CompareDistributions({x1,x4,x2,x5,x3},'tst','ranksum','lab',{'CON','RO','CIL','FENO','ISO'},'dsp','bar','ori',...
    'vert','tit',noms{ii},'sig','*', 'uni',uni{ii},'fig',0,...
    'fs',FS,'dc',[c11;c44;c22;c33;c55],'wi',0.5); axis square;
end
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOL2 '\trajectories_RyRs_cells.png']));



% MASS vs waves
FOLL='D:\Dades Lab\RyR-GFP_Activation_Sparks';
SS=dir([FOLL '\*RyR-GFP*']);
massSpkt=[]; SPKt=[];
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


% MASS vs FREQ
c33=[255, 204, 0]/255;
c55=[0, 134, 179]/255;
c11=[115, 0, 153]/255;
c44=[153, 0, 0]/255;
c22=[102, 153, 255]/255;

mg=0.1;
tit={'CON','RO','CIL','FENO','ISO'};
clf;
v=[1,2,3,4,5,11];
noms={'RyRs activated','RyRs initiators','Mean activations RyRs initiators',...
    'Max path length','max #activations per RyRs','PathLength·#RyRactivations'};
ii=1; v=[6];
indd={icon;iro;icil;ifeno;iiso};
cc={c11,c44,c22,c33,c55};
for jj=1:5
    subNM(2,3,jj,mg);
    x1=SPKt(indd{jj},v(ii)).*SPKt(indd{jj},21)+SPKt(indd{jj},v(ii)+6).*SPKt(indd{jj},22)...
        +SPKt(indd{jj},v(ii)+12).*SPKt(indd{jj},23);
    x1f=MMt(indd{jj},3);
    % x1f=FreqR(indd{jj},1);
    scatter(x1f,x1,50,cc{jj},'filled'); hold on;
    b1=x1f\x1; x=[0:max(x1f)/50:max(x1f)]; c=quantile(x1,0.05);
    plot(x1f,b1*x1f+c,':','Color',[0.5 0.5 0.5]);
    xlabel('<act RyR2 iniciators>'); ylabel('Mass (a.u.)');
    title(tit{jj});
end
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOLL '\traj_mass_actI.png']));


clf;
ii=1; v=[6];
indd={icon;iro;icil;ifeno;iiso};
for jj=1:5
    subNM(2,3,jj,mg);
    x1=SPKt(indd{jj},v(ii)).*SPKt(indd{jj},21)+SPKt(indd{jj},v(ii)+6).*SPKt(indd{jj},22)...
        +SPKt(indd{jj},v(ii)+12).*SPKt(indd{jj},23);
    x1f=MMt(indd{jj},4);
    % x1f=FreqR(indd{jj},1);
    scatter(x1f,x1,50,cc{jj},'filled'); hold on;
    b1=x1f\x1; x=[0:max(x1f)/50:max(x1f)]; c=quantile(x1,0.05);
    plot(x1f,b1*x1f+c,':','Color',[0.5 0.5 0.5]);
    xlabel('Path Length'); ylabel('Mass (a.u.)');
    title(tit{jj});
end
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOLL '\traj_mass_PL.png']));


clf;
ii=1; v=[6];
indd={icon;iro;icil;ifeno;iiso};
for jj=1:5
    subNM(2,3,jj,mg);
    x1=SPKt(indd{jj},22)...
        +SPKt(indd{jj},23);
    x1f=FreqR(indd{jj},1);
    scatter(x1f,x1,50,cc{jj},'filled'); hold on;
    b1=x1f\x1; x=[0:max(x1f)/50:max(x1f)]; c=quantile(x1,0.05);
    plot(x1f,b1*x1f+c,':','Color',[0.5 0.5 0.5]);
    xlabel('Spark Freq. (spk/(s·\mum^{2}))'); ylabel('fract >1RyRs spk');
    title(tit{jj});
end
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOLL '\traj_clustspk_spkfreq.png']));

clf;
ii=1; v=[6];
indd={icon;iro;icil;ifeno;iiso};
for jj=1:5
    subNM(2,3,jj,mg);
    x1=SPKt(indd{jj},22)...
        +SPKt(indd{jj},23);
    x1f=FreqR(indd{jj},2);
    scatter(x1f,x1,50,cc{jj},'filled'); hold on;
    b1=x1f\x1; x=[0:max(x1f)/50:max(x1f)]; c=quantile(x1,0.05);
    plot(x1f,b1*x1f+c,':','Color',[0.5 0.5 0.5]);
    xlabel('Waves Freq. (waves/s)'); ylabel('fract >1RyRs spk');
    title(tit{jj});
end
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOLL '\traj_clustspk_wavesfreq.png']));


clf;
ii=1; v=[20];
indd={icon;iro;icil;ifeno;iiso};
for jj=1:5
    subNM(2,3,jj,mg);
    x1=SPKt(indd{jj},v(ii)).*SPKt(indd{jj},21)+SPKt(indd{jj},v(ii)).*SPKt(indd{jj},22)...
        +SPKt(indd{jj},v(ii)).*SPKt(indd{jj},23);
    x1f=SPKt(indd{jj},6).*SPKt(indd{jj},21)+SPKt(indd{jj},6+6).*SPKt(indd{jj},22)...
        +SPKt(indd{jj},6+12).*SPKt(indd{jj},23);
    scatter(x1f,x1,50,cc{jj},'filled'); hold on;
    b1=x1f\x1; x=[0:max(x1f)/50:max(x1f)]; c=quantile(x1,0.05);
    plot(x1f,b1*x1f+c,':','Color',[0.5 0.5 0.5]);
    ylabel('<#clust/spk>');xlabel('Mass (a.u.)');
    title(tit{jj});
end
set(gcf,'position',[1 41 1920/2.5 964/2]);
saveWysiwyg(gcf,([FOLL '\traj_clustspk_mass.png']));
