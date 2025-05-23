function [trajectories]=EventTrajectories(EVfilt,folder,RF,Iryr,pos)

extra4='tajectories_spk';
if (exist([folder '\' RF '\' extra4],'dir')==0), mkdir([folder '\' RF '\' extra4]);end

trajectories={};
uu=unique([EVfilt{:,10}]); uu(uu==0)=[];
for ii=1:length(uu)
    ind=find([EVfilt{:,10}]==uu(ii)); % Troba events que comparteixen etiqueta
    if length(ind)>1,
        tt1=[EVfilt{ind,6}]; tt=tt1; % Ordre activació basat en el temps.... (MAL?)
        % Mirar si hi ha RyRs que s'activen al mateix moment
        [u,~,c]=unique(tt);
        if length(u)==1||max(c)==1, % Si s'activen a l'hora
            Apre=ones(length(ind),length(u));
            Apost=zeros(length(ind),length(u));
            C=Apost-Apre;
            c=zeros(length(tt),1);
        else
            dt=u(2:end)-u(1:end-1); % transicions
            Apre=zeros(length(ind),length(u)-1);
            Apost=zeros(length(ind),length(u)-1);
            Cpre=zeros(length(ind),length(u));
            vact=cell(length(dt),2);
            for jj=1:max(c)-1
                %ipre=unique(u(c(c==jj))); % Busca t1
                ipre=unique(u(c(c<=jj))); % Busca t1
                ipost=unique(u(c(c==jj+1))); % Busca t2
                indpre=[];
                for kk=1:length(ipre);
                    indprei=find(tt1==ipre(kk));indpre=[indpre, indprei];
                    % if jj==1, Cpre(kk,1)=1; end % Primers en activar-se
                end
                indpost=find(tt1==ipost);
                if length(indpre)==1, % en t1 només hi ha un activat --> activador
                    Apre(indpre,jj)=1;
                    if jj==1, Cpre(indpre,1)=1; end
                    for kk=1:length(indpost),
                        Apost(indpost(kk),jj)=1;
                        ind111=EVfilt{ind(indpre),1}; % RyR de indpre
                        ind222=EVfilt{ind(indpost(kk)),1}; % RyR de indpost
                        Cpre(indpost(kk),jj+1)=ind111;
                        vact{jj,1}=[vact{jj,1};[ind111,ind222]];
                        vact{jj,2}=[vact{jj,2};[[pos(ind111,1),pos(ind111,2);pos(ind222,1),pos(ind222,2)]]];
                    end
                else % Hi ha varios RyRs activats en t1
                    % Busca distancies:
                    x=[EVfilt{ind(indpre),4}]; y=[EVfilt{ind(indpre),5}];
                    for kk=1:length(indpost)
                        x2=[EVfilt{ind(indpost(kk)),4}]; y2=[EVfilt{ind(indpost(kk)),5}];
                        Z = squareform(pdist([x2,x;y2,y]','euclidean'));
                        zz=Z(1,:); zz(zz==0)=[]; [~,cc]=min(zz);
                        Apre(indpre(cc),jj)=1; % jj+1?
                        Apost(indpost(kk),jj)=1;
                        
                        %                         % REVISAR!!
                        %                         ipre=ind(indpre(cc)); ipost=ind(indpost(kk));
                        %                         v{jj}=[v{jj};[EVfilt{ipre,4},EVfilt{ipre,5}];[EVfilt{ipost,4},EVfilt{ipost,5}]];
                        %
                        % if jj==1, Apost(indpre(indpre~=cc),jj+1)=1;end
                        
                        ind111=EVfilt{ind(indpre(cc)),1};
                        ind222=EVfilt{ind(indpost(kk)),1};
                        Cpre(indpost(kk),jj+1)=ind111;
                        if jj==1, Cpre(indpre,1)=1;end;
                        vact{jj,1}=[vact{jj,1};[ind111,ind222]];
                        vact{jj,2}=[vact{jj,2};[[pos(ind111,1),pos(ind111,2);pos(ind222,1),pos(ind222,2)]]]; end
                    
                end
                
            end
            C=Apost-Apre; % C activation; tt dt; ind ROIs;
            trajectories{ii,4}=dt;
            trajectories{ii,5}=vact;
            trajectories{ii,6}=Cpre;
        end
        trajectories{ii,1}=C; trajectories{ii,2}=[tt',c,[EVfilt{ind,1}]'];trajectories{ii,3}=ind;  trajectories{ii,7}=[EVfilt{ind,1}];
    end
    trajectories{ii,3}=ind; % Events que paritcipen
    poss=[EVfilt{ind,4};EVfilt{ind,5}]';
    %     D=squareform(pdist(poss))*DX;
    
    
end

% Representar fletxes:
fact=2;
nf=4; % pix fletxa
cmapg=gray(256); cmapg(:,[1,3])=0; cmapg(end,:)=[1 1 1];
Iryr2=(Iryr-min(min(Iryr)))/(max(max(Iryr))-min(min(Iryr))); % Escalar Iryr de [0 : 1]
Iryr2=Iryr2.*.9; Iryr3=imresize(Iryr2,fact);
% figure(1),clf;set(gcf,'Position',[scz(1) scz(2) round(scz(3)/2) round(scz(4)/2)]);
for ii=1:length(uu)
    % cla;
    Iryr2=Iryr3;
    if ~isempty(trajectories{ii,1}), % ~isempty(trajectories{ii,5})
        try,
            if ~isempty(trajectories{ii,7})&&isempty(trajectories{ii,5}), % S'han coactivat al moment
                ryrs=trajectories{ii,7};
                for jj=1:length(ryrs) % Pinta els punts dels RyRs actius en el primer instant de temps
                    Iryr2(pos(ryrs(jj),1)*fact,pos(ryrs(jj),2)*fact)=1;
                end
                [allx,ally] = puntsEnmig([pos(ryrs(1),2)*fact pos(ryrs(2),2)*fact],[pos(ryrs(1),1)*fact pos(ryrs(2),1)*fact]);
                for xx=1:length(allx); Iryr2(ally(xx),allx(xx))=1; end;
                imwrite(Iryr2*256,cmapg,[folder,'\' RF '\' extra4 '\EventTraectories_' num2str(EVfilt{trajectories{ii,3}(1,1),10}) '.png']);
            end
            if ~isempty(trajectories{ii,5}),
                [a,b]=size(trajectories{ii,5}); % num transicions
                vact=trajectories{ii,5};
                ryrs=trajectories{ii,7};
                pre=trajectories{ii,6};
                indp=find(pre(:,1)>0);
                for jj=1:length(indp) % Pinta els punts dels RyRs actius en el primer instant de temps
                    Iryr2(pos(ryrs(indp(jj)),1)*fact,pos(ryrs(indp(jj)),2)*fact)=1;
                end
                Cul=[]; Cap=[];
                for jj=1:a
                    [c,d]=size(vact{jj,1}); % numero de parelles en una mateixa transició de temps
                    for kk=1:c
                        p2=vact{jj,2}(kk*2,:);p1=vact{jj,2}(kk*2-1,:);dp=p2-p1;
                        
                        [allx,ally] = puntsEnmig([p2(2)*fact p1(2)*fact],[p2(1)*fact p1(1)*fact]);
                        for xx=1:length(allx); Iryr2(ally(xx),allx(xx))=1; end;
                        
                        Cul=[Cul;[p1(2) p1(1)]*fact]; Cap=[Cap;[p2(2) p2(1)]*fact]; tamany=nf; col='w'; gruix=1;
                    end
                end
                
                if ~isempty(Cul),
                    for kk=1:size(Cul,1),
                        [ori1,des1,ori2,des2]=pintaCapCul(Cul(kk,:),Cap(kk,:),tamany,col,gruix);
                        ori1=round(ori1); ori2=round(ori2); des1=round(des1);des2=round(des2);
                        [allx,ally] = puntsEnmig([des1(1) ori1(1)],[des1(2) ori1(2)]);
                        for xx=1:length(allx); Iryr2(ally(xx),allx(xx))=1; end;
                        [allx,ally] = puntsEnmig([des2(1) ori2(1)],[des2(2) ori2(2)]);
                        for xx=1:length(allx); Iryr2(ally(xx),allx(xx))=1; end;
                    end;
                end
                % imagesc(Iryr2); axis image off; colormap(cmapg);% hold on;
                % title(['Event ' num2str(EVfilt{trajectories{ii,3}(1,1),24})]);
                %saveas(1,[folder,'\' RF '\' extra4 '\EventTraectories_' num2str(EVfilt{trajectories{ii,3}(1,1),24}) '.png']);
                %saveWysiwyg(1,[folder,'\' RF '\' extra4 '\EventTraectories_' num2str(EVfilt{trajectories{ii,3}(1,1),24}) '.png']);
                imwrite(Iryr2*256,cmapg,[folder,'\' RF '\' extra4 '\EventTraectories_' num2str(EVfilt{trajectories{ii,3}(1,1),10}) '.png']);
            end
        end
    end
end
close all;

end