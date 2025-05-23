function [events,signals,ev,eventsW]=EventsWavDetection(folder,RF,tww,Tww,win,signals,DT,rr,umbw,umbwW,nnr,pos,Iryr,V2,posr,scz,Ch,mam)
% Detecta events amb la wavelet i umbral

signalso=signals;

SIG=sum(signals)/rr;

noise=12; %ms
[infnoise,~]=escollirEscalesW(DT,noise,tww); % Per filtrar soroll
[inf,sup]=escollirEscalesW(DT,tww,Tww);

ev=zeros(size(signals));events=ev; eventsW=ev; eventsW2=ev;
wwct=0;
% Busca APs
s=zeros(1,round(nnr/10));s(1,round(nnr/10/2))=1;
ww=cwt(s,sup,'gaus2');
sap=zeros(rr,length(s));
for ii=1:rr, sap(ii,:)=ww; end
signalsAP=conv2(signals,sap,'same'); % imagesc(signalsAP);
eventsW3=double(signalsAP>39); % imagesc(eventsW3);

% if mean(mean(signals))<0.15&&mean(mean(signals))>0.09, umbwW=4; end % Cas RyR 31!
for kk=1:rr
    sig=signals(kk,:);
    % Noise
    wn=cwt(sig-mean(sig),infnoise,'gaus2'); sig(wn>1)=mean(sig); %   w=cwt(sig-mean(sig),[inf:sup],'gaus2');
    % plot(sig,'LineWidth',2,'Color','k');hold on; plot(wn(1,:),'LineWidth',2); % plot(w(end,:)); w=sum(w,1);
    
    % Events / Sparks
    w=cwt(sig-mean(sig),inf,'gaus2'); % plot(sig);hold on; plot(w);
    [MM,mm]=findPeaks6(w); % MM(w(MM)<mean(w(w>umbw))*.95)=[];
    MM(w(MM)<umbw)=[];
%            % plot(sig,'Color','k','LineWidth',2);
%             hold on;plot(w,'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]); plot(ones(1,nnr)*umbw,'Color',[0.9290, 0.6940, 0.1250]); %scatter(MM,w(MM));
%             % plot(ones(nnr,1)*umbw,'LineWidth',2);title('Sparks Detection');xlabel('Frames'); ylabel('Fnorm');
%             %legend('Fnorm','Wavelet filt','Sparks threshold');
%             plot(ones(round(1000/DT))*-.5,'Color','k','LineWidth',2);
%             text(3,-.3,['1000 ms'],'Color','k','FontSize',12); axis off;
%             set(gcf,'Position',[scz(1) scz(2) scz(3)/1.8 scz(4)/4]); %plot(signals(kk,:));
    
    % Waves
    ww=cwt(sig-mean(sig),sup,'gaus2');
    %     plot(sig,'LineWidth',2);hold on;plot(ww,'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]);
    %     plot(ones(nnr,1)*umbwW,'LineWidth',2,'Color',[0.4940, 0.1840, 0.5560]);title('Waves Detection');xlabel('Frames'); ylabel('Fnorm');
    %     legend('Fnorm','Wavelet filt','Waves threshold');  set(gcf,'Position',[scz(1) scz(2) scz(3)/2 scz(4)/4]);
    
    inz=find(ww>umbwW); % 1
    %     if ~isempty(inz),
    %         c=inz(1); v=find(diff(inz)~=1);int=inz(v); int=sort([int, inz(v+1)]); int=[c, int]; % Busca on s'inicia i on s'acaba la wave
    %         ind=[]; ind=find(ww>umbwW);
    %         if ~isempty(ind),
    %             wwct=wwct+1;
    %             if length(unique(diff(ind)))>2, % hi ha més d'una wave
    %                 cd=ind(1); v=find(diff(ind)>1); intd=ind(v+1); intd=[cd, intd];
    %             else
    %                 intd=min(ind);
    %             end
    %             for jj=1:length(intd)
    %                 v=int-intd(jj); v(v>0)=1; v(v<=0)=0;
    %                 in=find(v==0,1,'last'); fn=find(v==1,1,'first'); % in=int(in); fn=int(fn);
    %                 eventsW(kk,int(in):int(fn))=1;
    %             end
    %             eventsW(kk,ind)=1;
    %         end
    %     end
    v=[];
    if ~isempty(inz),
        wwct=wwct+1;
        v=find(diff(inz)>round(37/DT)); % Busca salts més grans que 5 ms
        if ~isempty(v),
            indd=[inz(v), inz(v+1)]; indd=sort(indd); indd=[inz(1), indd, inz(end)]; % Agafa el primer punt i els punts d'intercanvi
            for ii=1:length(indd)-1
                vv=indd(ii):indd(ii+1);
                if mean(ww(vv))>=umbwW,    eventsW(kk,vv)=1; end % El tram és més elevat en la wave
            end
        else
            eventsW(kk,inz)=1;
        end
    end
    
    eventsW2(kk,:)=ww;
    
    % Buscar punts màxims events/spk
    for jj=1:length(MM)
        mxw=MM(jj);
        ind1=max([1,mxw-round(win/2)]); ind2=min([length(sig),mxw+win]);
        ss1=sig(ind1:ind2);
        ind11=max([1,ind1-3]);ind22=max([ind2,length(sig)]);
        if all(ev(kk,ind11:ind22)==0),
            [~,realmax]=max(sig(ind1:ind2)); events(kk,realmax+ind1-1)=1;
            [~,realmin]=min(w(realmax+ind1-1:ind2)); ev(kk,realmax+ind1-1:realmin+realmax+ind1-1)=ones(length(realmax+ind1-1:realmin+realmax+ind1-1),1);
        end
    end
    
end
if size(ev,2)>nnr, ev=ev(1:end,1:nnr);end;

eventsWc=conv2(eventsW,[0 1 0;1 1 1;0 1 0],'same');
eventsW=double(eventsWc>0);
L=bwlabel(double(eventsWc>0)); % imagesc(L);
%  subplot(411);imagesc(signals); subplot(412);imagesc(eventsW2);
%  subplot(413);imagesc(eventsW3);subplot(414);imagesc(eventsW);linkaxes
if mam~=1,
    for ii=1:max(max(L))
        [r,c]=find(L==ii);
        if length(unique(r))<8, eventsW(L==ii)=0; end % si hi ha més de 6 RyRs activats, sinó macro spark (estava a 40)
    end
else
    for ii=1:max(max(L))
        [r,c]=find(L==ii);
         if length(unique(r))<8, eventsW(L==ii)=0; end % si hi ha més de 6 RyRs activats, sinó macro spark
    end
end
%  if wwct<8, eventsW(eventsW>0)=0;end % Si només hi ha 6 RyRs amb waves és falç --> macrosparks
eventsW=double((eventsW+eventsW3)>0); % imagesc(eventsW);
% subplot(211);imagesc(signals); subplot(212); imagesc(eventsW2>0.2);


save([folder,'\',RF,'\ev_det.mat'],'events','signals','pos','ev','eventsW');

% Representa:
rrr=rr;
representaPeaksWav(folder,RF,signals,events,ev,V2,pos,posr,Iryr,scz,nnr,rr,rrr,DT,Ch,eventsW,eventsW2)

end



% % PROVES:
%     % Buscar punts màxims events/spk
%     for jj=1:length(MM)
%         mxw=MM(jj);
%         ind1=max([1,mxw-round(win/2)]); ind2=min([length(sig),mxw+win]);
%         ss1=sig(ind1:ind2);
%         ind11=max([1,ind1-3]);ind22=max([ind2,length(sig)]);
%         if all(ev(kk,ind11:ind22)==0),
%             %                     ssd=ss1(2:end)-ss1(1:end-1);% ssd=[0,ssd];
%             %                     ssd(ssd<0)=0; [~,relmax]=max(ssd);
%             %                     zz=find(ssd==0); indz=zz-relmax; [~,z1]=max(indz(indz<0));[~,z2]=min(indz(indz>0));
%             %                     vv=ss1(zz(z1):zz(z1+z2)); [~,realmax]=max(vv);realmax=realmax+zz(z1)-1;
%             %                     %                     vv=min([length(ss1) realmax+1]):min([length(ss1) realmax+3]);
%             %                     %                     if length(vv)>1,if any(ssd(vv)>0.1), ind=find(ssd(vv)>0.1); realmax=realmax+ind(end); end,end
%             %                     events(kk,realmax+ind1-1)=1;
%             %                     %[~,realmax]=max(diff(ss1,2)); events(kk,realmax+ind1)=1;
%             [~,realmax]=max(sig(ind1:ind2)); events(kk,realmax+ind1-1)=1;
%             [~,realmin]=min(w(realmax+ind1-1:ind2)); ev(kk,realmax+ind1-1:realmin+realmax+ind1-1)=ones(length(realmax+ind1-1:realmin+realmax+ind1-1),1);
%         end
%     end