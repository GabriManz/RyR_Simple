function [] = RyRSimple(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% RyR Simple %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General settings
N=1;            % Number of experiments to process (set N=0 to analyze all sub-folders in folder)
NI=1;            % Number of images to process --> N==1!!!
tagRyR='ch00';  % relevant channel identifier
RF='ResRyR_prova_IMs_CNIC';    % name for results folder
sk=0;           % skip detection
zoom=0;         % Zoom image

% Detection settings
th=.01;        % noise factor
fr=1;           % number of frames to join (56)
tn=1;           % total number of resulting images
pg=2.5;           % Pealing rings (um)

% Filtering parameters
roiR=.5;          % roi radius around RyRs candidate center (in um)
thI=.05;           % intensity threshold
minR=.05;         % minimum radius (in um)
maxR=1.2;           % maximum radius (in um)
sr=1;             % show rejected sparks in output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% outputImage = RyRSimple(N,tagRyR,resultsf,skipdet,threshold,numberofframes,...
% numberofoutputimg,pealingrings,roiR,intensity,minR,maxR,showreject, folder);

RFGR=[];

if(nargin>0),
    N=varargin{1};
    NI=varargin{17};
    tagRyR=varargin{2};
    RF=varargin{3};
    sk=varargin{4};
    zoom=varargin{16};         % Zoom image
    
    th=varargin{5};
    fr=varargin{6};
    tn=varargin{7};
    pg=varargin{8};
    
    roiR=varargin{9};
    thI=varargin{10};
    minR=varargin{11};
    maxR=varargin{12};
    sr=varargin{13};
    
    RFGR=varargin{15}; % Ruta RyRGRcomp --> fer img pealing
    
    folder=varargin{14};
    [rutes,N]=inicialitza(N,folder);
else
    [rutes,N]=inicialitza(N);
end

%load('ultimpaquet.mat');



NN=N;
S=[];
warning('on','all');
if N==1&&NI>1,
    %N=NI;
    tn=NI;
    ensenya(['Starting with folder ' rutes{1}]);
    folder=rutes{1};
    [S,folder1] = uigetfile([folder '\*' tagRyR '.tif'],['Select ' num2str(NI) ' images'],'MultiSelect', 'on');
end

for ii=1:N
%     try
%         if NN==1,
            cc=1;
%         else
%             cc=ii;
%         end
        ensenya(['Starting with folder ' rutes{cc}]);
        % INICIALITZA
        folder=rutes{cc};

        if sk==0,
            ensenya(['Reading.']);
            [DX,DT,DZ] = buscaResolucio(folder,tagRyR);
            s=[dir([folder,'\*.tif']),dir([folder,'*.jpg']),dir([folder,'*.png'])];
            if length(s)<5, tn=1; end % Si només hi ha una imatge de dos canals (no serie de stacks)
            if length(s)<tn, tn=length(s); end % Si només hi ha una imatge de dos canals (no serie de stacks)
            
            %if NI==1,
                [Isum,Ch]=composeIms(folder,tagRyR,fr,tn,S);
%             else
%                 folder=folder1;
%                 [Isum,Ch]=composeIms(folder,tagRyR,fr,tn,S);
%             end
%             
            ensenya(['ROI size: ',num2str(1+2*ceil(roiR/DX)),'x',num2str(1+2*ceil(roiR/DX)),' pix.']);
            
            % DETECTA
            ensenya(['Processing.']);
            [IQsum, Q, Inorm]=processIms(Isum, DX, th,Ch);
            
            % mesura:Q=[X,Y,IntMitja,IntMax,R]
            ensenya(['Measuring.']);
            as=zlineangle(Q);
            % tic
            [Q, IZ, pealing]=mesures(Inorm, Q, ceil(roiR/DX), DX,as,IQsum, pg);
            % toc;
            if (exist([folder,'\',RF],'dir')==0), mkdir([folder,'\',RF]);end
            save([folder,'\',RF,'\detection.mat'],'Q','IZ','pealing','IQsum','Inorm', 'Isum', 'as','Ch','DX','DT','DZ');
        else
            %try
            load([folder,'\',RF,'\detection.mat']);
            %end
            %             if exist(Q)==0, % No s'ha fet el càlcul prèviament
            %                 sk=1;
            %                 RyRSimple(N,tagRyRg,RFg,sk,thg,fr,tn,pg,roiRg,thIg,minRg,maxRg,sr,folder,RFGR,zoom);
            %             end;
        end
        
        
        % FILTRA
        ensenya(['Filtering.']);
        [Qcom,Q,Qr]=filtraryrs(Q,Inorm, thI, minR, maxR); % diffR, minR i maxR en um
        % aveuredeteccio(Isum,Qcom,1);
        
        
        % GUARDA
        ensenya(['Saving.']);
        % Data
        guardaCSVs(Qcom, folder, RF, th, fr, tn, roiR, thI, minR, maxR, sr, as);
        % Images
        representaIm(Qcom, Inorm, Ch, DX, folder, RF,sr,zoom);
        %         if ~isempty(RFGR),
        %         if (exist([folder,'\',RF],'dir')==0), mkdir([folder,'\',RFGR]);end
        %         end
        RFGR=[];
        [ImP,cmap]=representaImP(Qcom,Q, Inorm, pealing, Ch, DX, folder, RF,RFGR,zoom);
        % clusterzlinesDFT(Q,DX,Inorm,folder,RF);
        
        %         try
        %             representaIm2(Qcom, Inorm, pealing, IZ, Ch, DX, folder, RF);
        %         end
        % representaZl(Qcom, Inorm, Ch, DX, folder, RF);
        
        % HTML
        ryrHtml(folder, RF, Q,Qr,th, fr, tn, roiR, thI, minR, maxR, as); % radi ROI en pix
        
        
%         if N>1,
%         Qcomi=Qcom; Qi=Q; Qri=Qr;
%         IZi=IZ; pealingi=pealing;
%         IQsumi=IQsum;Inormi=Inorm; Isumi=Isum;
%         ImPi=ImP; cmapi=cmap;
%         
%         Qcom={}; Q={}; Qr={};
%         IZ={}; pealing={};
%         IQsum={};Inorm={}; Isum={};
%         ImP={}; cmap={};
%         
%         Qcoms{ii}=Qcomi; Qs{ii}=Qi; Qrs{ii}=Qri;
%         IZs{ii}=IZi; pealings{ii}=pealingi;
%         IQsums{ii}=IQsumi;Inorms{ii}=Inormi; Isums{ii}=Isumi;
%         ImPs{ii}=ImPi; cmaps{ii}=cmapi;
%         
%         Qcom=Qcoms; Q=Qs; Qr=Qrs;
%         IZ=IZs; pealing=pealings;
%         IQsum=IQsums;Inorm=Inorms; Isum=Isums;
%         ImP=ImPs; cmap=cmaps;
%         end
        
        save([folder,'\',RF,'\finaldetection.mat'],'Qcom','Q','Qr','IZ','pealing','IQsum','Inorm','Isum','ImP','cmap','as','Ch','DX','DT','DZ');
        save([folder,'\',RF,'\parameters.mat'],'N','tagRyR','th','fr','tn','pg','roiR','thI','minR','maxR','sr');
        
%     catch me
%         warning(['Error, skipping experiment: ' me.message]);
%     end
    
    
end

end
