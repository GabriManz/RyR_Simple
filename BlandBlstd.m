function [blS,blSstd]=BlandBlstd(ev,eventsW,DT,rr,signals,nnr)
% Elimina waves i events de les senyals, estima el bl on no hi ha events.

% Bl and sigma:
EE=ev+eventsW; EE(EE>0)=1;
se = strel('line',round(250/DT),180);
EE2 = imdilate(EE,se,'same'); EE2(EE2>0)=1;
blS=zeros(1,rr); blSstd=blS;
for ii=1:rr
    evv=EE2(ii,:);
    sig=signals(ii,:); sig(evv>0)=[];
    blS(ii)=quantile(sig,.1); % Base line
    blSstd(ii)=1.4826*mad(sig); % noise?
%     plot(signals(ii,:)); hold on; plot(ones(length(signals(ii,:)))*(blS(ii)+2.5*blSstd(ii)));
    % plot(signals(ii,:)); hold on; plot(ones(length(signals(ii,:)))*(blS(ii)+4*blSstd(ii)));
%     pause;
%     hold off;
end
        
end
        