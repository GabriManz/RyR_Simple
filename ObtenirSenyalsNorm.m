function [pos,posr,Q,Iryr,rr,rrr,sign, signals,scz,lim,signalsM,limmad]=ObtenirSenyalsNorm(G,folder, RF,DX,nnr,rs,a,b,V,DT,rp)
% Obten senyals dels RyRs
% Normalitza senyals respecte el bl local

% RyRs i noRyRs
Q=G.Q{1};
pos=[Q(:,3),Q(:,2)]; rrr=length(pos); % y=pos(ii,1); (files) x=pos(ii,2); (columnes)
posr=pos;
Iryr=G.Isum;
%[posnor]=trobarPuntsbtwZl(pos,G,DX,rrr,folder,RF); close all; pos=[posr;posnor];

%%%%%%%%%%%%%%%%%% OBTENIR SENYALS
rr=length(pos); sign=zeros([length(pos),nnr]);
r=ceil(rs/DX);
% ROIs quadrades:
for ii=1:rr
    for jj=1:nnr
        mny=max([1 pos(ii,1)-r]); mxy=min([pos(ii,1)+r a]);
        mnx=max([1 pos(ii,2)-r]); mxx=min([pos(ii,2)+r b]);
        sign(ii,jj)=mean(mean(V(mny:mxy,mnx:mxx,jj)));
    end
end


%         % Escala/normalitza % (F-Fo)/Fo:
%         signals=zeros(size(sign));
%         lim=quantile2(sign(:),.01); mm=zeros(rr,1);
%
%         if lim*1.5<mean(mean(sign)), lim=quantile2(sign(:),.05); end;
%         % kk=65; kk=82;me=mean(mean(sign)); mt=quantile(sign(:),001);
%         for kk=1:rr
%             F=sign(kk,:); %     if mean(F)<me*.7,  F=F+me*.5; end
%             % Fo=max([lim,quantile2(F,.01)]); % if Fo==0, Fo=0.0001; end % quantile2(F,.01)
%             Fo=lim;
%             sig=((F-Fo)/Fo);
%             signals(kk,:)=sig;
%         end
%         % Save signals:
%         save([folder,'\',RF,'\signals.mat'],'sign','signals','pos','posr');

% Signals:
scz = get( 0, 'Screensize' );
figure(1); clf;
for kk=1:rr
    plot(sign(kk,:)); hold on;
end
xlabel('Frames'); ylabel('F'); ylim([0 80]);
title(['F mean ' num2str(quantile(sign(:),.5)),' max = ',num2str(max(sign(:))),' q5 = ',num2str(quantile(sign(:),.05))]);
set(gcf,'Position',[scz(1),scz(2),scz(3)/3,(scz(4))/4]);
saveWysiwyg(1,[folder '\' RF '\Sign_all.png']);

% Escala/normalitza % (F-Fo)/mad(Fo):
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
    SIGm = medfilt1(SIG,ww,'truncate'); % plot(SIGm);
    M = movstd(SIG,ww); % plot(SIGm); hold on; plot(M); plot(SIGm+M*3)
    Mm=mean(M);
    if (Mm*2)+quantile(F,.05)>mean(F), lim3=median(SIGm+M*3);SIG2=SIG;SIG2(SIGm+M*3>lim3)=[]; 
    else
    % Anterior 190829:
    lim3=median(SIGm+M*3);SIG2=SIG;SIG2(SIG2>lim3)=[];
    end
%     % New 190829:
%     ind=find(M<quantile(M,0.15)); % plot(F(ind)-M(ind))
%     SIG2=F(ind);
    Fbasal=[Fbasal,SIG2];
    
%     q05=quantile(SIGm,.05);
%     q95=quantile(SIGm,.95);
%     
%     %Fo=median(SIGm-M*2); lim(kk)=Fo;
%     if(Mm*2)+quantile(F,.05)>mean(F), Fo=median(SIG2);
%         if ((q95-q05)/3)<(q95-Fo), Fo=quantile(SIG2,0.15); end
%     else, Fo=median(SIGm);
%     end
    Fo=median(SIGm); 
    lim(kk)=Fo; Fo2=Fo; % plot(SIGm); hold on; plot(ones(1,nnr)*Fo);
    limmad(kk)=mad(SIG2); %sig=((F-Fo)/Fo2);
    % if Fo>20, Fo=25; end; % Motiu QQ plot! %% 181122
    if FF<=4, Fo2=Fo*(5-Fo); end % Arreglar Baseline més petit que 1!
    sig=((F-Fo)/Fo2)*10;
    signals(kk,:)=sig;
%     sig=((F-Fo))*10;
%     signals(kk,:)=sig;
    AF(kk,:)=(F-Fo);
    signalsM(kk,:)=(F-Fo)/mad(SIG2);
end


% 191016 RyRs no existents
[r,c]=find(isinf(signals));
rr2=unique(r);
pos(rr2,:)=[]; posr=pos;
signals(rr2,:)=[]; sign(rr2,:)=[];
signalsM(rr2,:)=[];
AF(rr2,:)=[];
limmad(rr2,:)=[]; lim(rr2,:)=[];
rr=length(pos); rrr=rr;


% limmad; mean(limmad);
signals=signals/30;
signalsM=signalsM/30;
pll=0;
if pll==1,
    scz = get( 0, 'Screensize' );
    plot(sign(kk,:),'LineWidth',3); hold on; plot(SIGm,'LineWidth',3); plot(M,'LineWidth',3); plot(ones(nnr,1)*Fo,'LineWidth',3);
    xlabel('Frames'); ylabel('F'); legend('F','Fmed','Fstd','Fo');set(gcf,'Position',[scz(1),scz(2),scz(3)/2,(scz(4))/4]);
    hold off;
    
    subplot(211), plot(sign(kk,:),'LineWidth',3); hold on;plot(ones(nnr,1)*lim3,'LineWidth',3); title('F');
    hold off;
    subplot(212), plot(SIG2,'LineWidth',3);xlabel('Frames');ylabel('F'); title('Fbasal');
    set(gcf,'Position',[scz(1),scz(2),scz(3)/2,(scz(4))/4]);
    hold off;
    
    subplot(211), plot(sign(kk,:),'LineWidth',3);xlabel('Frames');ylabel('F'); title('F');
    subplot(212), plot(signals(kk,:),'LineWidth',3);ylim([0 1]);xlabel('Frames');ylabel('F'); title('Fnorm');
    set(gcf,'Position',[scz(1),scz(2),scz(3)/2,(scz(4))/4]);
end
% plot(sum(signals)/rr);
% scz = get( 0, 'Screensize' );
% plot(signals(kk,:)); xlabel('Frames');ylabel('(F-Fo)/MAD(Fo)'); set(gcf,'Position',[scz(1),scz(2),scz(3)/3,(scz(4))]);

figure(1), clf;
subplot(211);imagesc(sign,[0 50]); title(['mean signals ' num2str(mean(mean(sign)))  ' max ' num2str(max(max(sign)))]);
subplot(212);imagesc(signals,[0 1]); title(['mean signalsN ' num2str(mean(mean(signals))) ' max ' num2str(max(max(signals)))]);
saveWysiwyg(1,[folder '\' RF '\provaNorm.png']);

figure(1), clf;
subplot(211);imagesc(sign,[0 50]); title(['mean signals ' num2str(mean(mean(sign)))  ' max ' num2str(max(max(sign)))]);
subplot(212);imagesc(signalsM,[0 1]); title(['mean signalsN ' num2str(mean(mean(signalsM))) ' max ' num2str(max(max(signalsM)))]);
saveWysiwyg(1,[folder '\' RF '\provaNormMAD.png']);

figure(1),clf;
[counts,centers] =hist(Fbasal,[0:max(max(Fbasal))/30:max(max(Fbasal))]);
bar(centers,100*counts/length(Fbasal),'FaceColor',[0.9290, 0.6940, 0.1250]); axis square; ylabel('%'); xlabel('Values Fbasal');
title(['F baseline ' num2str(mean(lim)),' mean MAD = ' num2str(mean(limmad))]); %set(gcf,'Position',[scz(1),scz(2),scz(3)/2,(scz(4))]);
saveWysiwyg(1,[folder '\' RF '\Fbasal_distrib_all.png']);

% Save signals:
save([folder,'\',RF,'\signals.mat'],'sign','signals','pos','posr','lim','Fbasal','limmad','signalsM');
%save([folder,'\',RF,'\signals.mat'],'sign','signals','pos','posr','lim','limmad','Fbasal');

% REPRESENTA NORM:
ind=find(folder=='\');
name=folder(ind(end)+1:end);
% Signals:
scz = get( 0, 'Screensize' );
figure(1); clf;
for kk=1:rr
    plot(AF(kk,:)); hold on;
end
xlabel('Frames'); ylabel('F'); ylim([-10 60]);
title(['Fo raw = ' num2str(mean(lim)),'AF mean = ' num2str(quantile(AF(:),.5)),' max = '...
    ,num2str(max(AF(:))),' q5 = ',num2str(quantile(AF(:),.05)) ' q95 = ',num2str(quantile(AF(:),.95))]);
set(gcf,'Position',[scz(1),scz(2),scz(3)/3,(scz(4))/4]);
saveWysiwyg(1,[folder '\' RF '\Signals_F-Fo_' name '.png']);

scz = get( 0, 'Screensize' );
figure(1); clf;
for kk=1:rr
    plot(signals(kk,:)); hold on;
end
xlabel('Frames'); ylabel('F'); ylim([0 3]);
title(['Fnorm mean ' num2str(quantile(signals(:),.5)),' max = '...
    ,num2str(max(signals(:))),' q5 = ',num2str(quantile(signals(:),.05)) ' q95 = ',num2str(quantile(signals(:),.95))]);
set(gcf,'Position',[scz(1),scz(2),scz(3)/3,(scz(4))/4]);
saveWysiwyg(1,[folder '\' RF '\Signals_Norm_all.png']);


%%%%%%%%%%%%%%%%%% REPRESENTA SENYALS:
scz = get( 0, 'Screensize' );
if rp==1,
    figure(1);clf,imagesc(signals); axis image; set(gcf,'Position',[scz(1),scz(2),scz(3),(scz(4)/2)]);
    title('RyR signals'); xlabel('Frames'); ylabel('F/Fo');saveWysiwyg(1,[folder,'\',RF,'\signals.png']);
    % sig=sum(signals,1); % plot(sig);
    
    delta=3.5;
    figure(2),clf, set(gcf,'Position',scz);
    for ii=1:rr, s=signals(ii,:)+ii*delta; plot(s); hold on; mm=mean(s);if ii<rrr,text(-30,mm,num2str(Q(ii,1))); else text(-30,mm,'No');end, end
    set(gca,'ytick',[]); axis([-50 2500 0 max(s)+10]); title('RyR signals'); xlabel('Frames'); ylabel('F/Fo');
    saveWysiwyg(2,[folder,'\',RF,'\RyR_signals.png']);
    close all;
end


end