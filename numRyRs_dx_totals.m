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
vsf=[1:16]; ctff=1; masst=[]; Dxt=[]; RR=[]; NZ=[]; Dit=[]; dnncellst=[];
for sss=1:length(vsf)%size(SS,1)
    FOLc=[FOLL '\' SS(vsf(sss)).name]; Sc=dir([FOLc '\*RyR-GFP*']);
    ind=find([Sc(:).isdir]==0);
    for ii=1:length(ind), Sc(ind(ii)-ii+1)=[]; end
    massSpk=zeros(size(Sc,1),4); % mass
    massSpk2=zeros(size(Sc,1),4); % mass
    SPK=zeros(size(Sc,1),26); % Sparks parameters
    nRyRs=zeros(size(Sc,1),4); % num RyRs / Events
    nRyRs2=zeros(size(Sc,1),12); % num RyRs / Events
    Dc=zeros(size(Sc,1),5); rrc=zeros(size(Sc,1),5); NZc=zeros(size(Sc,1),6);
    dnncells=zeros(size(Sc,1),5);
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
            
             Di=[]; rri=[]; rri2=[]; NZi=[]; dnnryrst=zeros(size(S,1),5);
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
                
                % Dist DNN RyRs
                dnnryrs=mean(G.Q{1,1}(:,7));
                dnnryrst(ss,1)=dnnryrs;
                dnnryrst(ss,2)=std(G.Q{1,1}(:,7))/sqrt(length(G.Q{1,1}(:,7)));
                dnnryrst(ss,3)=mean(G.Q{1,1}(:,6));
                dnnryrst(ss,4)=std(G.Q{1,1}(:,6))/sqrt(length(G.Q{1,1}(:,6)));
                dnnryrst(ss,5)=indx;
                % dnnryrst=[dnnryrst;dnnryrst];
                
                freqi=0; densi=0; 
                if ~isempty(EVfilt),
                    signalso=signals;
                    
                    et=unique([EVfilt{:,10}]); % Sparks
                    pos2=zeros(size(distcov,1),4); D=[];
                    for ii=1:size(distcov,1),
                        pos2(ii,1:2)=[distcov{ii,2}];
                        indr=find([EVfilt{:,10}]==et(ii));
                        ind=length(unique([EVfilt{indr,1}]));
                        indt=[EVfilt{indr,6}];
                        pos2(ii,3)=min(indt);
                        pos2(ii,4)=ind; % num RyRs
                        try,
                            D(ii,1)=distcov{ii,1}(1);
                            D(ii,2)=distcov{ii,1}(2);
                            D(ii,3)=ind;
                        end
                    end
                    Di=[Di;D];
                    % distcov=[D cent CovM [Sx,Sy] eigV sqrt(lambda) ve*lamb Iev];
                    %          1  2     3     4      5       6          7     8
                    % angRyRs=[SZ DZ nRyrs nSZ nDZ v]
                    %           1  2   3    4   5  6
                    % aRyRs=[nRyRs tRyRs SZ DZ v amp RoR t2p DHM tau]
                    %          1     2   3  4  5  6   7   8   9   10
                    ind=find(angRyRs(:,3)==2);
                    if ~isempty(ind),
                        dx=zeros(length(ind),2);
                        for ii=1:length(ind), 
                            [Ms]=aRyRs{ind(ii),5};
                            ind1=find(Ms(:,2)==1);
                            ind2=find(Ms(:,2)==2);
                            if ~isempty(ind1), dx(ii,1)=min([distcov{ind(ii),1}]); 
                                if dx(ii,1)>1, dx(ii,1)=dx(ii,1)/1.5; end; if dx(ii,1)>1,dx(ii,1)=0.8;end 
                            end
                            if ~isempty(ind2), dx(ii,2)=max([distcov{ind(ii),1}]);
                                if dx(ii,2)<1.6, dx(ii,2)=dx(ii,2)*1.2; end; % if dx(ii,2)<1.6, dx(ii,2)=1.7; end;
                                % if dx(ii,2)>2.2, dx(ii,2)=2.06; end;
                            end
                        end
                        dx1=0; dx2=0;
                        try, dx1=dx(:,1); dx1(dx1==0)=[]; dx1=mean(dx1); end
                        try, dx2=dx(:,2); dx2(dx2==0)=[]; dx2=mean(dx2); end
                        dx1(isnan(dx1))=0;  dx2(isnan(dx2))=0;
                            NZi=[NZi;sum(angRyRs(ind,4)), sum(angRyRs(ind,5)),size(EVfilt,1),dx1,dx2];
                    end
                    
                    ind=length(unique([EVfilt{:,1}]));
                   
                    rri2=[rri2;ind/rr];
                end
                 rri=[rri;rr];
            end
            dnncells(ssc,1:5)=mean(dnnryrst,1);
            
            if size(Di,1)>1,
            Di(Di(:,3)~=2,:)=[];
            Dc(ssc,1)=mean(Di(:,1)); Dc(ssc,2)=mean(Di(:,2)); Dc(ssc,3)=mean(Di(:,3));
            Dc(ssc,4)=indx;
            Dc(ssc,5)=sss;
            end
            Dit=[Dit;[Di,ones(size(Di,1),1)*indx]];
            
            rrc(ssc,1)=round(mean(rri)); 
            if ~isempty(rri2),rrc(ssc,2)=round(round(mean(rri))*max(rri2));
            rrc(ssc,3)=mean(rri2); end
            rrc(ssc,4)=indx;
            rrc(ssc,5)=sss;
            
            if ~isempty(NZi), 
                NZc(ssc,1)=sum(NZi(:,1)); 
                NZc(ssc,2)=sum(NZi(:,2));
                NZc(ssc,3)=sum(NZi(:,3));
                NZc(ssc,4)=indx; 
                ind=find(NZi(:,4)>0);
                if ~isempty(ind),NZc(ssc,5)=mean(NZi(ind,4)); end % dx SZ
                ind=find(NZi(:,5)>0);
                if ~isempty(ind), NZc(ssc,6)=mean(NZi(ind,5)); end % dx DZ
            end
        
        end
    end
    save([FOLc '\dx.mat'],'Dc','rrc','NZc')
    Dxt=[Dxt;Dc];
    RR=[RR;rrc]; % nº RyRs i nº RyRs actius
    NZ=[NZ;NZc];
    dnncellst=[dnncellst;dnncells];
end
save([FOLL '\dx_ryrs.mat'],'Dxt','RR','NZ','Dit')
save([FOLL '\dnn_ryrs.mat'],'dnncellst');

icon=find(NZ(:,4)==1);
ifeno=find(NZ(:,4)==5);
iiso=find(NZ(:,4)==3);
Tcon = array2table(NZ(icon,[1:3,5:6]),'VariableNames',{'nSZ','nDZ','ntotals','dxSZ','dxDZ'});
Tfeno = array2table(NZ(ifeno,[1:3,5:6]),'VariableNames',{'nSZ','nDZ','ntotals','dxSZ','dxDZ'});
Tiso = array2table(NZ(iiso,[1:3,5:6]),'VariableNames',{'nSZ','nDZ','ntotals','dxSZ','dxDZ'});
writetable(Tcon,[FOLL '\dx_2cluste_ev_CON.xlsx']);
writetable(Tfeno,[FOLL '\dx_2cluste_ev_FENO.xlsx']);
writetable(Tiso,[FOLL '\dx_2cluste_ev_ISO.xlsx']);
% Tcon = array2table(Dxt(icon,:),'VariableNames',{'dx','cent','nRyRs','drug','mouse'});
% Tfeno = array2table(Dxt(ifeno,:),'VariableNames',{'dx','cent','nRyRs','drug','mouse'});
% Tiso = array2table(Dxt(iiso,:),'VariableNames',{'dx','cent','nRyRs','drug','mouse'});
% writetable(Tcon,[FOLL '\dnn_ryrs_CON.xlsx']);
% writetable(Tfeno,[FOLL '\dnn_ryrs_FENO.xlsx']);
% writetable(Tiso,[FOLL '\dnn_ryrs_ISO.xlsx']);


icon=find(dnncellst(:,5)==1);
ifeno=find(dnncellst(:,5)==5);
iiso=find(dnncellst(:,5)==3);
Tcon = array2table(dnncellst(icon,:),'VariableNames',{'DNN','SE_dnn','radi','SE_radi','drug'});
Tfeno = array2table(dnncellst(ifeno,:),'VariableNames',{'DNN','SE_dnn','radi','SE_radi','drug'});
Tiso = array2table(dnncellst(iiso,:),'VariableNames',{'DNN','SE_dnn','radi','SE_radi','drug'});
writetable(Tcon,[FOLL '\dnn_ryrs_CON.xlsx']);
writetable(Tfeno,[FOLL '\dnn_ryrs_FENO.xlsx']);
writetable(Tiso,[FOLL '\dnn_ryrs_ISO.xlsx']);

Dxt(isnan(Dxt(:,1)),:)=[];
icon=find(Dxt(:,4)==1);
ifeno=find(Dxt(:,4)==5);
iiso=find(Dxt(:,4)==3);
iro=find(Dxt(:,4)==4);
icil=find(Dxt(:,4)==2);

RyRdist=0.9;
clf;
icond=find(Dit(:,4)==1);
[h,c]=hist(Dit(:,1),[0:0.05:2.5]);
counts=[0:0.05:2.5]; counts2=[0:0.005:2.5]
s1=h/sum(h);  s1=interp1(counts,s1,counts2);
windowSize = 20;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
f1=filtfilt(b,a,s1); f1=f1/max(f1); f1=f1*max(h/sum(h));
plot(counts2,f1,'LineWidth',3,'Color',[143, 172, 199]/255);
hold on;
y=[0:0.05:0.15];
plot(ones(length(y),1)*RyRdist,y,':','Color',[0.6,0.6,0.6],'LineWidth',3);
set(gca,'box','off','fontsize',16);
xlim([0 3]);
saveWysiwyg(gcf,(['D:\Dades Lab\RyR-GFP_Activation_Sparks\FIGURAS\200223\dist_RyRs_spks.png']));




% Distància entre RyRs de sparks de 2 RyRs
dx=mean(Dxt(icon,1)); % dist entre sparks de 2 RyRs
dx=mean(Dxt(iro,1)); % dist entre sparks de 2 RyRs
dx=mean(Dxt(icil,1)); % dist entre sparks de 2 RyRs
dx=mean(Dxt(ifeno,1)); % dist entre sparks de 2 RyRs
dx=mean(Dxt(iiso,1)); % dist entre sparks de 2 RyRs



% nombre de RyRs i nº RyRs activats
icon=find(RR(:,4)==1);
ifeno=find(RR(:,4)==5);
nRyRst=sum(RR(icon,1));
nRyRsta=sum(RR(icon,2));

nRyRst=sum(RR(ifeno,1));
nRyRsta=sum(RR(ifeno,2));

% Salts línes Z
icon=find(NZ(:,4)==1);
iro=find(NZ(:,4)==4);
icil=find(NZ(:,4)==2);
ifeno=find(NZ(:,4)==5);
iiso=find(NZ(:,4)==3);

SZ=sum(NZ(iiso,1))
DZ=sum(NZ(iiso,2))
SZ/(SZ+DZ)
DZ/(SZ+DZ)

% NZ = [nSZ nDZ nspk drug dSZ dDZ];
%        1   2    3   4    5   6
icon=find(NZ(:,4)==1);
ifeno=find(NZ(:,4)==5);
iiso=find(NZ(:,4)==3);
iro=find(NZ(:,4)==4);
icil=find(NZ(:,4)==2);

ind={}
ind{1}=icon;
ind{2}=icil;
ind{3}=iro;
ind{4}=ifeno;
ind{5}=iiso;
TT=zeros(5,4);
for ii=1:5
    % dist
    dSZi=NZ(ind{ii},5); dSZi(dSZi==0)=[];
    dDZi=NZ(ind{ii},6); dDZi(dDZi==0)=[];
    dSZ=mean(dSZi)
    dDZ=mean(dDZi)
    SEdSZ=std(dSZi)/sqrt(length(dSZi))
    SEdDZ=std(dDZi)/sqrt(length(dDZi))
    TT(ii,1)=dSZ; TT(ii,2)=SEdSZ;
    TT(ii,3)=dDZ; TT(ii,4)=SEdDZ;
  
end

for jj=1:size(NZ,1)
    NZ(jj,1:2)=NZ(jj,1:2)/max(sum(NZ(jj,1:2)));
    
end

TT2=zeros(5,4);
for ii=1:5
    
       % num salts
     dSZi=NZ(ind{ii},1); 
    dDZi=NZ(ind{ii},2); 
    
    dSZ=mean(dSZi)
    dDZ=mean(dDZi)
    SEdSZ=std(dSZi)/sqrt(length(dSZi))
    SEdDZ=std(dDZi)/sqrt(length(dDZi))
     TT2(ii,1)=dSZ; TT2(ii,2)=SEdSZ;
    TT2(ii,3)=dDZ; TT2(ii,4)=SEdDZ;
  
    
end

CONT=NZ(ind{1},6); CONT(CONT==0)=[];
for ii=2:5
    DRUGt=NZ(ind{ii},6);DRUGt(DRUGt==0)=[];
    [p,h]=ranksum(CONT,DRUGt)
    TT3(ii)=p;
end

CONT=NZ(ind{1},1); %CONT(CONT==0)=[];
for ii=2:5
    DRUGt=NZ(ind{ii},1);%DRUGt(DRUGt==0)=[];
    [p,h]=ranksum(CONT,DRUGt)
    TT3(ii)=p;
end


% SZ
for ii=1:5
    for jj=1:5
    DRUGt1=NZ(ind{ii},1);%DRUGt(DRUGt==0)=[];
    DRUGt2=NZ(ind{jj},1);
    [p,h]=ranksum(DRUGt1,DRUGt2)
    TT4(ii,jj)=p;
    end
end

% DZ
for ii=1:5
    for jj=1:5
    DRUGt1=NZ(ind{ii},2);%DRUGt(DRUGt==0)=[];
    DRUGt2=NZ(ind{jj},2);
    [p,h]=ranksum(DRUGt1,DRUGt2)
    TT5(ii,jj)=p;
    end
end