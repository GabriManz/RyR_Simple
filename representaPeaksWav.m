function  representaPeaksWav(folder,RF,signals,events,ev,V2,pos,posr,Iryr,scz,nnr,rr,rrr,DT,Ch,eventsW,eventsW2)

ind=find(folder=='\');

figure(1); clf;
subplot(311); imagesc(signals);title('RyRs signals'); xlabel('Frames'); ylabel('F/Fo');
subplot(312); imagesc(events); title('Peaks activation'); xlabel('Frames'); ylabel('Detection');
subplot(313); imagesc(ev); title('Activation'); xlabel('Frames'); ylabel('Detection');
linkaxes;
set(gcf,'Position',scz);
saveWysiwyg(1,[folder,'\',RF,'\Detection.png']);

figure(1); clf;
subplot(211); imagesc(signals);title('RyRs signals'); xlabel('Frames'); ylabel('F/Fo');
subplot(212); imagesc(events);title('Peaks activation'); xlabel('Frames'); ylabel('Detection');
linkaxes;
set(gcf,'Position',scz);
saveWysiwyg(1,[folder,'\',RF,'\Peaks.png']);
% saveWysiwyg(1,[folder,'\',RF,'\Zoom_Saltos.png']);

figure(1);clf;
subplot(311); imagesc(signals);title('RyRs signals'); xlabel('Frames'); ylabel('F/Fo');
subplot(312); imagesc(eventsW2, [0 5]);title('wavelet'); xlabel('Frames'); ylabel('Detection');
subplot(313); imagesc(eventsW);title('waves activation'); xlabel('Frames'); ylabel('Detection');
linkaxes;
set(gcf,'Position',scz);
saveWysiwyg(1,[folder,'\',RF,'\Waves_' folder(ind(end)+1:end) '.png']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Representa vídeos %%%%%%%%%%%%%


V3=V2.*0.6;
V4=V3;
% VV=zeros(size(V3));
for ii=1:nnr
    for jj=1:rr
        if ev(jj,ii)==1,
            V3(pos(jj,1),pos(jj,2),ii)=1;
            % VV(pos(jj,1),pos(jj,2),ii)=1;
            V4(pos(jj,1),pos(jj,2),ii)=1;
        else
            if jj<rrr,
                V4(posr(jj,1),posr(jj,2),ii)=1-5*(1/256);
            end
        end
    end
end
cm=gray(256); cm(1:154,[2:3])=0; cm(155:255,[1,3])=0; cm(end,:)=[1 1 1];
ind=find(folder=='\');
namefile=['video_act_RyRs_',folder(ind(end)+1:end)];
makemovie6(V3,[folder '\' RF],DT,1,1,5,0, Ch, namefile, 0);
namefile=['video_act_RyRs&GFP_',folder(ind(end)+1:end)];
makemovie62(V4,[folder '\' RF],DT,1,1,5,0, Ch, namefile, 0,cm); % Representa RyRs


% for ii=1:nnr, V4(:,:,ii)=ind2rgb(V4(:,:,ii),cm);end
% V4=zeros([size(V3,1),size(V3,2),3]);
% V4(:,:,:,1)=V3;
% for ii=1:nnr,V4(:,:,ii,2)=G.IQsum;end;
% V4(:,:,:,3)=VV;
% V4=permute(V4,[1 2 4 3]);
% clf, figure(1),for ii=1:nnr, imagesc(V3(:,:,ii),[0 1]); colormap(cm); axis image; pause(0.1); end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fer representacions 3D %%%%%%%%
% - RyRs and singals
% - RyRs and activation


% % Signals
% figure(88),clf;
% for ii=1:rr
%     x=ones(1,nnr)*pos(ii,2);
%     y=10*signals(ii,:)+pos(ii,1);
%     t=1:1:nnr;
%     plot3(t,x,y);
%     hold on;
% end
% [a, b]=size(Iryr);
% cmap=gray(256); cmap(:,[1,3])=0;
% xImage = [0 0; 0 0];   % The x data for the image corners
% yImage = [0 b; 0 b];   % The y data for the image corners
% zImage = [0 0; a a];   % The z data for the image corners
% surf(xImage,yImage,zImage,'CData',Iryr,'FaceColor','texturemap'); colormap(cmap);
% xlabel('Frames');
% zlabel('F/Fo');
% title('RyRs and signals');
% view([20 -20 20]);
% saveWysiwyg(88,[folder,'\',RF,'\RyRsandSignals3D.png']);
% savefig(88,[folder '\' RF '\RyRsSignals.fig']);
% 
% 
% % % Activation (És molt lent l'scatter)
% % figure(86),clf;
% % for ii=1:rr
% %     [xx,yy]=find(events(ii,:)==1);
% %     if ~isempty(xx),
% %         for jj=1:length(xx)
% %             x=pos(ii,2);
% %             y=pos(ii,1);
% %             t=yy(jj);
% %             if ii<rrr, col='.g'; else, col='.r'; end;
% %             scatter3(t,x,y,40,col);
% %             hold on;
% %         end
% %     end
% % end
% % [a, b]=size(Iryr);
% % cmap=gray(256); cmap(:,[1,3])=0;
% % xImage = [0 0; 0 0];   % The x data for the image corners
% % yImage = [0 b; 0 b];   % The y data for the image corners
% % zImage = [0 0; a a];   % The z data for the image corners
% % surf(xImage,yImage,zImage,'CData',Iryr,'FaceColor','texturemap'); colormap(cmap);
% % xlabel('Frames');
% % zlabel('F/Fo');
% % title('Activation');
% % view([20 -20 20]);
% % saveWysiwyg(86,[folder,'\',RF,'\RyRsandDet3D.png']);

end