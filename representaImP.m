function [Im,cmap]=representaImP(Qcom,Q, Inorm, pealing, Ch, DX, folder, RF,RFGR,zoom)

fact=1;
if zoom==0, fact=ceil(.5/DX); end 


ind=find(folder=='\',1,'last');
file=folder(ind+1:end);

for ii=1:length(Qcom)
    filename=['RyRsPealing',num2str(ii),'_',file];
    %I=zeros(size(Inorm(:,:,1)));
    Qi=Qcom{ii};
    Qii=Q{ii};
    
    IQ=zeros(size(Inorm));
    for jj=1:length(Qii)
        IQ(Qii(jj,3),Qii(jj,2))=1;
    end
    
    Im=Inorm(:,:,ii); 
    mx2=max(max(pealing(:,:,ii)));
%     mx=0;
%     mx2=0;
%     for kk=2:length(Qi)
%         % if(mx<Qi{kk,7}),mx=Qi{kk,7};end
%         if(mx2<Qi{kk,8}),mx2=Qi{kk,8};end % pealing
%     end

    % Normalitza img segons el colormap que se li assignarà
    top=255-mx2;
    %     top=256-(mx+1+1)-2;
    % cont=0;
    Im=(top-2)*(Im-min(min(Im)))/(max(max(Im))-min(min(Im)));
    % Im=enxufaLlegenda2(Im,DX,1); %%% ARREGLAR!!
    Im(Im>top-1)=top;
    Im(IQ>0)=top;
    
    Im=imresize(Im,fact,'nearest');
    

    % Colormap pealing
    % cmp=gray(mx2+2); cmp(2,:)=[];
    cmp=hsv(mx2);
    
    % Colormap cell
    cmapc=gray(top+1);
    vect=[1:size(cmapc,2)]; vect(vect==Ch)=[]; cmapc(:,vect)=0; % Es queca amb el canal del Ch
    cmapc(end,:)=[1 1 1];
    
    % Unió colormaps
    %     cmap=[cmap;cmn;cmz;cmp];
    cmap=[cmapc;cmp];
    
    
    
   
    
    %     % Posa num de RyR
    %     for kk=2:length(Qi)
    %         % Posar nums
    %         color=top+1;
    %         if sr==1,
    %             if ~strcmp(Qi{kk,9},'Ok'),
    %                 color=top;
    %             end
    %             Im=textIm(Qi{kk,1}*fact,Qi{kk,2}*fact,num2str(kk-1),Im,'textcolor',color,'blending','off');
    %         else
    %             if(strcmp(Qi{kk,9},'Ok')),
    %                 Im=textIm(Qi{kk,1}*fact,Qi{kk,2}*fact,num2str(kk-1),Im,'textcolor',color,'blending','off');
    %             end
    %         end
    %     end
    
    
    % Posar pealing
    pealing2=imresize(pealing(:,:,ii),fact,'nearest');
    Ip=zeros(size(pealing2)); % Ip2=Ip;
    % romans={'I' 'II' 'III' 'IV' 'V' 'VI' 'VII' 'VIII' 'IX' 'X'};
    for kk=1:max(max(pealing2))
        color=kk+top;
        aquests=pealing2>=kk; aquests=double(aquests);%imagesc(aquests);
        I2=imdilate(aquests,[0 1 0; 1 1 1; 0 1 0]);
        resta=I2-aquests; resta(resta==1)=color;
        Ip=Ip+resta;
        [a,b]=find(resta(1,:)>0);
        if isempty(b), % Si no son zooms de cell
            [~,b]=find(resta(round(end/2),:)>0);
            a=round(size(resta,1)/2);
            %             if isempty(b), % Si no son zooms de cell
            %                  [a,b]=find(resta>0);
            %                  b(1)=b(round(length(bb)/2));
            %                  a(1)=a(round(length(bb)/2));
            %             end
        end
        try Ip=textIm(b(1),a(1)+10,[' ' romanNumerals(kk)],Ip,'textcolor',color,'horizontalalignment','left','blending','off');
%             if kk~=max(max(pealing2)),
%                 Ip=textIm(b(2),a(1)+10,[romanNumerals(kk) ' '],Ip,'textcolor',color,'horizontalalignment','right','blending','off');
%             end;
        end
    end
    % imagesc(Ip);
    % Im(Ip>0)=256;
    [ind]=find(Ip>0);
    for kk=1:length(ind)
        Im(ind(kk))=Ip(ind(kk));
    end
    
%     for kk=1:mx
%         color=kk+top+1;
%         Im=textIm(bt(kk)-2,at(kk),[num2str(kk)],Im,'textcolor',color,'horizontalalignment','right','blending','off');
%     end
    % imagesc(Im);colormap(cmap); axis image;
    
    
    % Im=imresize(Im,fact,'nearest');
    
    
    imwrite(uint8(Im),cmap,[folder,'\',RF,'\',filename,'.png']);
%     if ~isempty(RFGR),
%     imwrite(uint8(Im),cmap,[folder,'\',RFGR,'\',filename,'.png']);
%     end;
    
end


end
