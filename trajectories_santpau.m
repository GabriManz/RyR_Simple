% TRAJECTORIES & Dist RyRs
close all;
clear all;
clc;

%%%%%%%%%%%%%%%%% FOLDER MOUSE EXPERIMENTS
ensenya('Select the experiment folder');
FOLc='D:\Dades Lab\RyR-GFP_prova_Sant_Pau\170215_RyR-GFP30'; % MOUSE FOLDER
% Folders normal strure:
%   > MICE (MVM_1) -->  select this one!
%       > CELLS (MVM_1_cell01)
%           > SERIES (MVM_1_cell01_serie1)

%%%%%%%%%%%%%%%%% Parameters
usuari=1;
factampt=6;
% Event grouping (SAME you had fot activation parameters)
wt1=350;                % Temporal window (ms) wt=20; 181004 wt = 50;
ccorr=0.5;              % Correlation
dspk1=2.5;%1;           % um

RF='ResRyR_12_2'; % Results folder



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sc=dir([FOLc '\*RyR-GFP*']);

SPKin=zeros(size(Sc,1),7);
SPKnoin=zeros(size(Sc,1),7);
SPKin2=zeros(size(Sc,1),7);
SPKnoin2=zeros(size(Sc,1),7);
RFo=RF;
factamp=factampt(1);  ctf=1; FOLtt={};
% spk in waves
SPKint=[]; SPKnoint=[]; SPKint2=[]; SPKnoint2=[];
for ssc=1:size(Sc,1) % Cell
    MMt=[];
%     if usuari==1,
%         FOL=uigetdir([FOLc],'Select the main experiments folder'); S=dir([FOL '\*RyR-GFP*']);
%     else
%         FOL=[FOLc '\' Sc(ssc).name]; S=dir([FOL '\*RyR-GFP*']);
%     end
    FOL=[FOLc '\' Sc(ssc).name]; S=dir([FOL '\*RyR-GFP*']);
    S([S.isdir]==0)=[];
    
    %%%% PREGUNTAR VERO
%     ind1=[];
%     for ss=1:size(S,1),
%         indls=strfind(S(ss).name,'LS');
%         imgse=strfind(S(ss).name,'Image');
%         sestim=strfind(S(ss).name,'Serie1');
%         if any([~isempty(indls),~isempty(imgse),~isempty(sestim)])
%             ind1=[ind1;ss];
%         end
%     end
%     if ~isempty(ind1), S(ind1)=[];end
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
        FOL=[FOLc '\' Sc(ssc).name]; % S=dir([FOL '\*RyR-GFP*']);
        
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
        for ii=1:size(Att2,1),Ra(ii)=Att2(ii,ii);end
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
        if ~isempty(PLt),mNP=max(min([max(Att)', PLt']'));else mNP=0; end
        
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
    
    % save
    T = array2table(MM,'VariableNames', {'RyRs_with_act','RyRs_acts_percent',...
        'Activacions_RyR_act','Max_path_length','max_node','max_PL_or_nodes','time_waves',...
        'nRyRs','drug','max_Path_Length_plus_percent_Activated','mean_Path_Length_plus_percent_Activated'});
    writetable(T,[FOL '\Trajectories_summary_cell.xlsx']);
end

% spk in wave
icon=find(SPKin(:,7)==1);
ifeno=find(SPKin(:,7)==5);
iiso=find(SPKin(:,7)==3);
Tcon = array2table(SPKin(icon,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
Tfeno = array2table(SPKin(ifeno,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
Tiso = array2table(SPKin(iiso,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
writetable(Tcon,[FOLc '\spk_inicio_waves_CON.xlsx']);
writetable(Tfeno,[FOLc '\spk_inicio_waves_FENO.xlsx']);
writetable(Tiso,[FOLc '\spk_inicio_waves_ISO.xlsx']);
icon=find(SPKnoin(:,7)==1);
ifeno=find(SPKnoin(:,7)==5);
iiso=find(SPKnoin(:,7)==3);
Tcon = array2table(SPKnoin(icon,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
Tfeno = array2table(SPKnoin(ifeno,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
Tiso = array2table(SPKnoin(iiso,:),'VariableNames',{'amp','RoR','t2p','FDHM','tau','nRyRs','drug'});
writetable(Tcon,[FOLc '\spk_no_inicio_waves_CON.xlsx']);
writetable(Tfeno,[FOLc '\spk_no_inicio_waves_FENO.xlsx']);
writetable(Tiso,[FOLc '\spk_no_inicio_waves_ISO.xlsx']);

