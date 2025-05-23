function [Vdet,labb]=RepresentaNoWavesClust(folder,RF,S,ss,factamp,DX,DT,G,nnr,scz,signals,V2,Ch,EVfilt,rr,delta,wst)

% CLUSTERS AND ACTIVATIONS
Qr=G.Q; Qr=Qr{1};
ind=unique([EVfilt{:,1}]);
labb=unique([EVfilt{:,10}]);
cmclust=superjet(length(labb),'lines');
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
saveWysiwyg(1,[folder,'\',RF,'\Signals_clust_AD_noWn_' S(ss).name '_' num2str(factamp) '.png']);

% Make Detection Volume NO WAVES:
Vdet=zeros(size(V2)); Vdete=Vdet; V3=V2*0.9; % V5=V2*(255-et);
% cmp=superjet(et,'Lines'); cmm=gray(256-et); cmm(:,2:3)=0; cm=[cmm;cmp];
fact=2;V4=zeros(size(V3,1)*fact, size(V3,2)*fact, size(V3,3)); ett=0;
if fact~=1, for ii=1:nnr, V4(:,:,ii)=imresize(V3(:,:,ii),fact,'nearest');end; end
for ii=1:size(EVfilt,1)
    x=EVfilt{ii,4}; y=EVfilt{ii,5}; ttsi=EVfilt{ii,6}; tts=[ttsi:min([ttsi+20 nnr])]; % !!
    
    Vdet(y,x,ttsi)=1;
    V3(y,x,tts)=1;
    Vdete(y,x,tts)=EVfilt{ii,10};% V5(y,x,tts)=(255-et)+EVfilt{ii,10};
    
    for jj=1:length(tts)
        I=imresize(V3(:,:,tts(jj)),fact,'nearest');
        if any(unique(V4(:,:,tts(jj)))~=0), I2=V4(:,:,tts(jj)); I2(I==1)=1; I=I2; end
        if ett~=[EVfilt{ii,10}], % ja s'ha posat l'etiqueta abans
            I=textIm(2+x*fact,2+y*fact,num2str([EVfilt{ii,10}]),I,'textcolor',[1 1 1],'blending','off');
        end
        V4(:,:,tts(jj))=I;
    end
    ett=[EVfilt{ii,10}];
end % figure; for ii=1:nnr, imagesc(V3(:,:,ii), [0 1]); axis image; pause(0.1); end
ind=find(folder=='\');namefile=['clust_actRyRs_noW',folder(ind(end)+1:end) '_Filt'];
makemovie6(V4,[folder '\' RF],DT,1,1,5,0, Ch, namefile, 0);


end