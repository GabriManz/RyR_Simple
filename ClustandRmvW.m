function [EVfilt,Vdet]=ClustandRmvW(folder,RF,S,ss,factamp,DX,DT,pos,G,wt,dspk,rmw,EVcorr,EV,nnr,rr,Iryr,events,eventsW,signals,ampF,tauF,rorF,t2pF,dhmF,r2F,V2,scz,blS,blSstd,Ch,limmad,lim)
% Hierarchical clust
% Uprise
% Representation
% Remove waves
% Representation without waves

ind=[];if size(EVcorr,1)>0, ind=find([EVcorr{:,19}]==1);  EV2=zeros(length(ind),6);end% indz=find([EVcorr{:,19}]==0);
EVfilt={}; EVfilt2=EVfilt; spkF=0;  Vdet=zeros(size(V2));
if ~isempty(ind),
    % Regenera les matrius ev info i ordena les etiquetes:
    for ii=1:length(ind),
        for jj=1:size(EVcorr,2),
            if jj<=20, EVfilt{ii,jj}=EVcorr{ind(ii),jj};
            else  EVfilt{ii,jj+6}=EVcorr{ind(ii),jj};
            end
        end,
    end
    for ii=1:length(ind),for jj=1:size(EV,2),EV2(ii,jj)=EV(ind(ii),jj);end,end;
    EVcorr2=EVcorr; EV=EV2;
    
    ensenya('Clustering sparks.');
    % [EVfilt,et,Iryr]=eventsRyRsClust3(wt,DT,DX,dspk,signals,Iryr,rr,nnr,EV,EVfilt);
    [EVfilt,et,Iryr]=eventsRyRsClust4(wt,DT,DX,dspk,signals,Iryr,rr,nnr,EV,EVfilt); % Hierarchical clust
    close all;
    % representaEventsiRyRs2(folder,RF,EVfilt, DT,Iryr,et,scz,nnr,signals,V2,pos,DX,rr)
    
    % Upraise calcium:
    [EVfilt]=UpRiseCalcium(EVfilt, DT,signals);
    save([folder,'\',RF,'\ClustEv_' num2str(factamp) '.mat'],'EVfilt','EVcorr','wt','DT','dspk','DX','events','signals','V2','pos','EV','blS');
    
    % Wave / no Wave:
    for ii=1:size(EVfilt,1),
        roi=EVfilt{ii,1}; punt=EVfilt{ii,6}; 
        if eventsW(roi,punt)>0, 
            EVfilt{ii,25}=1; 
        else
            EVfilt{ii,25}=0; 
        end; 
    end
    
    % FILTERED AND ACCEPTED
    delta=.25;
    wst=round(65/DT);
    RepresentaFiltered(EVfilt,DX,DT,delta,wst,rr,signals,events,eventsW,EVcorr,factamp,folder,RF,scz,nnr,blS,blSstd,S,ss,V2,limmad,lim);
    % representaEventsiRyRs2(folder,RF,EVfilt, DT,Iryr,et,scz,nnr,signals,V2,pos,DX,rr)
    
    % Make Detection Volume:
    Vdet=zeros(size(V2)); Vdete=Vdet; V3=V2*0.9; % V5=V2*(255-et);
    V32=V2.*0.6; % cmp=superjet(et,'Lines'); cmm=gray(256-et); cmm(:,2:3)=0; cm=[cmm;cmp];
    fact=2;V4=zeros(size(V3,1)*fact, size(V3,2)*fact, size(V3,3)); V42=V4;
    for ii=1:nnr,for jj=1:rr,V32(pos(jj,1),pos(jj,2),ii)=1-5*(1/256);end;end % GFP
    if fact~=1, for ii=1:nnr, V4(:,:,ii)=imresize(V3(:,:,ii),fact,'nearest'); V42(:,:,ii)=imresize(V32(:,:,ii),fact,'nearest');end; end
    for ii=1:size(EVfilt,1)
        x=EVfilt{ii,4}; y=EVfilt{ii,5}; tts=EVfilt{ii,6}; tts=[tts:min([tts+15 nnr])]; % !!
        Vdet(y,x,tts)=1; V3(y,x,tts)=1; V32(y,x,tts)=1;
        Vdete(y,x,tts)=EVfilt{ii,10};% V5(y,x,tts)=(255-et)+EVfilt{ii,10};
        for jj=1:length(tts)
            I=imresize(V3(:,:,tts(jj)),fact,'nearest');
            I=textIm(2+x*fact,2+y*fact,num2str([EVfilt{ii,10}]),I,'textcolor',[1 1 1],'blending','off');
            V4(:,:,tts(jj))=I;
            % GFP
            I2=imresize(V32(:,:,tts(jj)),fact,'nearest');
            I2=textIm(2+x*fact,2+y*fact,num2str([EVfilt{ii,10}]),I2,'textcolor',[1 1 1],'blending','off');
            V42(:,:,tts(jj))=I2;
        end
    end % figure; for ii=1:nnr, imagesc(V3(:,:,ii), [0 1]); axis image; pause(0.1); end
    ind=find(folder=='\');namefile=['clust_actRyRs_',folder(ind(end)+1:end) '_Filt'];
    makemovie6(V4,[folder '\' RF],DT,1,1,5,0, Ch, namefile, 0);
    namefile=['clust_actRyRs_GFP_',folder(ind(end)+1:end) '_Filt'];
    cm=gray(256); cm(1:154,[2:3])=0; cm(155:255,[1,3])=0; cm(end,:)=[1 1 1];
    makemovie62(V42,[folder '\' RF],DT,1,1,5,0, Ch, namefile, 0,cm); % Representa RyRs
    
    
    V3=V2.*0.6;
    V4=V3;
    % VV=zeros(size(V3));
    for ii=1:nnr
        for jj=1:rr
            %if ev(jj,ii)==1, % commented 190307
                V3(pos(jj,1),pos(jj,2),ii)=1;
                % VV(pos(jj,1),pos(jj,2),ii)=1;
%                 V4(pos(jj,1),pos(jj,2),ii)=1;
            % else
%                 if jj<rrr,
                    V4(pos(jj,1),pos(jj,2),ii)=1-5*(1/256);
%                 end
            % end
        end
    end
    cm=gray(256); cm(1:154,[2:3])=0; cm(155:255,[1,3])=0; cm(end,:)=[1 1 1];
    ind=find(folder=='\');
    namefile=['video_act_RyRs_',folder(ind(end)+1:end)];
    makemovie6(V3,[folder '\' RF],DT,1,1,5,0, Ch, namefile, 0);
    namefile=['video_act_RyRs&GFP_',folder(ind(end)+1:end)];
    makemovie62(V4,[folder '\' RF],DT,1,1,5,0, Ch, namefile, 0,cm); % Representa RyRs
    
    
    % TREURE WAVES
    EVfilt2=EVfilt;
    if rmw==1,
        [EVfilt2,EVfilt]=removeWaves(EVfilt,EVfilt2);
    end
    
    if ~isempty(EVfilt), % tots els sparks estaven en una wave
        [Vdet,labb]=RepresentaNoWavesClust(folder,RF,S,ss,factamp,DX,DT,G,nnr,scz,signals,V2,Ch,EVfilt,rr,delta,wst);
        
        % Spark frequency:
        sf=zeros(1,nnr); % et=length(unique([EVfilt{:,10}]));
        for ii=1:length(unique([EVfilt{:,10}])), ind=find([EVfilt{:,10}]==labb(ii)); sf(EVfilt{ind(1),6})=1; end
        spkF=numel(sf==1)/(rr*nnr*DT); % spk/(ms*RyRs)
        
        % EVcorr: [ROI ROIryr R/nR x y t sig(wt) evind indevcorr et amp AMP ror t2p dhm tau r2 bline in fin];
        % ind:      1     2     3  4 5 6   7      8        9     10  11  12  13  14 15  16  17   18   21 22
        % fextra=[]; sparks=[];
        save([folder,'\',RF,'\RyRsandEv_' num2str(factamp) '.mat'],'events','Vdet','signals','EVcorr','EVfilt','EVfilt2','spkF','ampF','tauF','rorF','t2pF','dhmF','r2F','V3');
        close all;
    else
        labb=[]; Vdet=V2;
        save([folder,'\',RF,'\RyRsandEv_' num2str(factamp) '.mat'],'events','Vdet','signals','EVcorr','EVfilt','EVfilt2','spkF','ampF','tauF','rorF','t2pF','dhmF','r2F','V3');
    end
else
    Vdet=zeros(size(V2));
    save([folder,'\',RF,'\RyRsandEv_' num2str(factamp) '.mat'],'events','Vdet','signals','EVcorr','EVfilt','EVfilt2','ampF','tauF','rorF','t2pF','dhmF','r2F');
end


end