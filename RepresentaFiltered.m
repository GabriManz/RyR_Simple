function []=RepresentaFiltered(EVfilt,DX,DT,delta,wst,rr,signals,events,eventsW,EVcorr,factamp,folder,RF,scz,nnr,blS,blSstd,S,ss,V2,limmad,lim)

% FILTERED AND ACCEPTED
tts=1:nnr; rr2=1:rr;
figure(1),clf;set(gcf,'Position',scz);

    % Wave / no Wave:
    for ii=1:size(EVfilt,1),roi=EVfilt{ii,1}; punt=EVfilt{ii,6}; if eventsW(roi,punt)>0, EVfilt{ii,25}=1; else  EVfilt{ii,25}=0; end; end
    
for jj=1:rr
    s=signals(jj,:)+jj*delta;
    plot(tts*DT,s,'k'); hold on;
    % plot([1:nnr]*DT,ones(nnr,1)'*(blS(jj)+4*blSstd(jj))+jj*delta);
    plot([1:nnr]*DT,ones(nnr,1)'*((limmad(jj)/lim(jj)*(1/3))*factamp*1.2)+jj*delta);
    inde=find(events(jj,:)==1);
    try, indef=find([EVfilt{:,1}]==jj);
        indew=find([EVcorr{:,1}]==jj);end
    if ~isempty(inde),
        for kk=1:length(inde)
            tt=[max([1 inde(kk)-wst]):min([nnr inde(kk)+wst+4])];
            
            % Det wavelet
            plot(tt*DT,s(tt),'Color','y');
            mm=mean(s);
            text((tts(1)*DT)-20,mm,num2str(rr2(jj)),'Color','k');
        end
        if ~isempty(indew),
            col2='r';
            for kk=1:length(indew)
                tt=[max([1 EVcorr{indew(kk),6}-15]):min([nnr EVcorr{indew(kk),6}+20])];
                if strfind([EVcorr{indew(kk),20}],'Low amp'), col2='r'; end
                if strfind([EVcorr{indew(kk),20}],'tau'), col2='m'; end
                if strfind([EVcorr{indew(kk),20}],'t2p'), col2='c'; end
                if strfind([EVcorr{indew(kk),20}],'border'), col2='g'; end
                if strfind([EVcorr{indew(kk),20}],'fdhm'),col2='c'; end
                
                % Det wavelet
                plot(tt*DT,s(tt),'Color',col2);
                mm=mean(s);
                % text((tts(1)*DT)-20,mm,num2str(rr2(jj)),'Color','k');
            end
        end
        try, for kk=1:length(indef)
                tt=[max([1 EVfilt{indef(kk),6}-15]):min([nnr EVfilt{indef(kk),6}+20])];
                
                % Det wavelet
                plot(tt*DT,s(tt),'Color','b'); end;
        end
    end
end
xlabel('t (ms)'); set(gca,'ytick',[]);
saveWysiwyg(1,[folder,'\',RF,'\Signals_events_' S(ss).name '_' num2str(factamp) '.png']);


% CLUSTERS AND ACTIVATIONS
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
            tt=[max([1 c-5]):min([nnr c+wst+5])];
            % Cluster
            plot(tt*DT,s(tt),'Color',cmclust(lab,:));
            
            % Act / Diff
            cool='.r'; if [EVfilt{inde(kk),25}]==1, cool='.y'; end % ETIQUETA wave
            scatter(c*DT,s(c)+.2,cool);
            scatter(EVfilt{inde(kk),23}*DT,s(EVfilt{inde(kk),23})+.2,15,'k','filled');
            %pause;
            mm=mean(s);
            text((tts(1)*DT)-20,mm,num2str(EVfilt{inde(kk),2}),'Color','k');
        end
    end
end
xlabel('t (ms)'); set(gca,'ytick',[]);
saveWysiwyg(1,[folder,'\',RF,'\Signals_clust_AD' S(ss).name '_' num2str(factamp) '.png']);

end