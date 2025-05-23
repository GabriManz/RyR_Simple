function [EVcorr,et,Iryr]=eventsRyRsClust4(wt,DT,DX,dspk,signals,Iryr,rr,nnr,EV,EVcorr)
% Fa matriu events amb tots els paràmetres (EVcorr)
% Clusteritza els events Hierarchical clistering
% Representa el resultats

wto=wt/DT; wt2=wto; % wt2=round(wto/2);
dspk2=dspk/DX;

% EV = [ ROI RyRroi ryr x  y  t ]
%         1     2    3  4  5  6

P = dspk2/wt2;
EV2=EV(:,4:6);
%EV2(:,1:2)=EV2(:,1:2)*DX; EV2(:,3)=EV2(:,3)*DT;
EV2(:,3)=EV2(:,3).*P; % Escalem la variable de temps

% DD=zeros(size(EV2,1));
% for ii=1:size(EV2,1)
%     for jj=1:size(EV2,1)
%         DD(ii,jj)=sqrt(((EV2(ii,1)-EV2(jj,1))^2)+((EV2(ii,2)-EV2(jj,2))^2)+((EV2(ii,3)-EV2(jj,3))^2));
%     end
% end


if size(EV,1)>1,
    % Z =linkage(EV2,'ward','euclidean'); % dendrogram(Z,'ColorThreshold',dspk2);
    D=pdist(EV2);
    Z =linkage(D);
    C = cluster(Z,'cutoff',dspk2,'Criterion','distance'); % C = cluster(Z,'cutoff',dspk2*dspk2,'Criterion','distance');
    %EV2(:,3)=EV2(:,3)./P;
    %scatter3(EV2(:,3),EV2(:,1),EV2(:,2),50,C,'filled');
    % for ii=1:size(EV,1),EVcorr{ii,10}=C(ii);end;
    
    % Ordena etiquetes del 1 al lab per ordre temporal:
    C2=C;
    for ii=1:length(C), ind=find(C==C(ii)); C2(ind)=ii; end
    u=unique(C2);
    for ii=1:length(u),ind=find(C2==u(ii));for jj=1:length(ind), C3(ind(jj))=ii; end,end
    C=C3;
    
    et=max(C);
    wt2=round(wt2);
    for ii=1:size(EV,1)
        tt=EV(ii,6); in=max([1 tt-wt2]); fin=min([tt+wt2 nnr]);
        s1=signals(EV(ii,1),in:fin); % Senyal en window temporal de l'event
        
        EVcorr{ii,7}=s1;
        EVcorr{ii,10}=C(ii);
    end
else
    tt=EV(1,6); in=max([1 tt-wt2]); fin=min([tt+wt2 nnr]);
    s1=signals(EV(1,1),in:fin);
    et=1;
    EVcorr{1,7}=s1;
    EVcorr{1,10}=1;
end



end