function [distcov, eryr]=distancesCovariance(EVfilt,Iryr,DX)

uu=unique([EVfilt{:,10}]); uu(uu==0)=[];
distcov=cell(length(uu),11);
eryr=zeros([8, length(uu)]);

for ii=1:length(uu)
    ind=find([EVfilt{:,10}]==uu(ii)); % Troba events que comparteixen etiqueta
    % if EVfilt{ind(1),3}==1,eryr(1,ii)=1;eryr(4,ii)=1;end; % S'inicia en un RyR (temp, no RoR)
    indd=ind;
    if length(indd)==1, % Nomes hi ha un RyR
        eryr(3,ii)=indd; % ROIs
        distcov{ii,2}=[EVfilt{indd,4},EVfilt{indd,5}];
        eryr(8,ii)=1;
        eryr(4,ii)=1;
    end
    if length((unique([EVfilt{indd,1}])))>1,
        % try
        Iev=Iryr;
        eryr(2,ii)=1; % Hi ha més d'un RyR activat
        eryr(4,ii)=length(indd); % Number of activated RyR
        % [~,c]=min([EVfilt{ind(indd),6}]); % temps %[~,c]=max([EVfilt{ind(indd),13}]); % RoR
        ror=[EVfilt{indd,13}]; ror=(ror-min(ror))/(max(ror)-min(ror));
        tt=[EVfilt{indd,6}]; eryr(5,ii)=max(tt)-min(tt); % dt peaks
        tt=(tt-min(tt))/(max(tt)-min(tt));tt=1-tt;
        first=tt+ror; [~,c]=max(first);
        eryr(3,ii)=indd(c); % First act
        % Variance of FW
        x=[EVfilt{indd,4}];y=[EVfilt{indd,5}];
        D=squareform(pdist([x;y]','euclidean'))*DX; D=D(:); D(D==0)=[];%sqrt(((x1-x2)^2)+((y1-y2)^2));
        eryr(6,ii)=mean(D);eryr(7,ii)=max(D);
        eryr(8,ii)=length(unique([EVfilt{indd,1}])); % número de RyRs
        
        c=[mean(x),mean(y)]; % Centroide
        C=cov(x,y); Sx=sqrt(C(1,1))*DX; Sy=sqrt(C(2,2))*DX; % Covariance matrix + variança x + variança y
        distcov{ii,1}=D; distcov{ii,2}=c; distcov{ii,3}=C; distcov{ii,4}=[Sx,Sy];
        [v,lamb]=eig(C); lam=unique(lamb); lam(lam==0)=[]; lam=sqrt(lam);
        distcov{ii,5}=v; distcov{ii,6}=lam; % sqrt(lambda)
        distcov{ii,7}=[Sx,Sy];
        
        % Eixos variança PCA:
        if length(indd)>2, % Si hi ha més de dos RyRs
            if length(lam)==1, lam(1)=lam; lam(2)=lam; end
            P1=[v(1,1)*lam(1),v(1,2)*lam(1)]; P2=[v(2,1)*lam(2),v(2,2)*lam(2)];
            distcov{ii,7}=[v(1,1)*lam(1),v(1,2)*lam(1);v(2,1)*lam(2),v(2,2)*lam(2)]*DX;
            P11=[-v(1,1)*lam(1),-v(1,2)*lam(1)]; P22=[-v(2,1)*lam(2),-v(2,2)*lam(2)];
            P1=round(P1+c); P11=round(P11+c); P2=round(P2+c); P22=round(P22+c);
            [xx1,yy1] = puntsEnmig([P1(1),P11(1)],[P1(2),P11(2)]);
            [xx2,yy2] = puntsEnmig([P2(1),P22(1)],[P2(2),P22(2)]);
            xx1(xx1<1)=[]; yy1(yy1<1)=[]; xx2(xx2<1)=[]; yy2(yy2<1)=[];
            xx1(xx1>size(Iev,2))=[];yy1(yy1>size(Iev,1))=[];
            xx2(xx2>size(Iev,2))=[]; yy2(yy2>size(Iev,1))=[];
            % Pinta:
            nn=min([length(xx1),length(yy1)]);
            nn2=min([length(xx2),length(yy2)]);
            nn3=min([length(x),length(y)]);
            for kk=1:nn,Iev(yy1(kk),xx1(kk))=max(max(Iryr))+3;end
            for kk=1:nn2,Iev(yy2(kk),xx2(kk))=max(max(Iryr))+3;end
            for kk=1:nn3,Iev(y(kk),x(kk))=max(max(Iryr))+1; end % Pinta RyRs event
            Iev(round(c(2)),round(c(1)))=max(max(Iryr))+2; 
           % Iev=enxufaLlegenda2(Iev,DX,max(Iev),'SE'); % Pinta centroide % imagesc(Iev);axis image;
            distcov{ii,8}=Iev;
            % Arregla cov:
            Sxi=zeros(1,length(indd)-1); Syi=zeros(1,length(indd)-1);
            for jj=1:length(indd)-1
                x=[EVfilt{indd(jj:jj+1),4}];y=[EVfilt{indd(jj:jj+1),5}];
                 C=cov(x,y);
                 Sxi(jj)=sqrt(C(1,1))*DX;
                 Syi(jj)=sqrt(C(2,2))*DX;
            end
            distcov{ii,4}=[max(Sxi),max(Syi)];
        end
        
        % Width:
        distcov{ii,12}=length((unique([EVfilt{indd,1}]))); % nº RyRs
        if ~isempty(distcov{ii,7}),
            if ~isempty(distcov{ii,8}), % Més de 2 RyRs
                P1=distcov{ii,7}(1,:);  P2=distcov{ii,7}(2,:); P11=-P1; P22=-P2;
                width1=sqrt(((P1(1)-P11(1))^2)+(P1(2)-P11(2))^2);
                width2=sqrt(((P2(1)-P22(1))^2)+(P2(2)-P22(2))^2);
                distcov{ii,10}=width1; distcov{ii,11}=width2; % spark gran
            else % Si són dos RyRs
                distcov{ii,10}=min(distcov{ii,7}); distcov{ii,11}=max(distcov{ii,7}); % x i y
            end
        else
            distcov{ii,10}=0.33; distcov{ii,11}=0.33; % Spark d'un RyR
        end
        
        %end
        clear v lamb lam C c x y Sx Sy P1 P2 P11 P22 xx1 xx2 yy1 yy2;
    end;
     distcov{ii,9}=length(unique([EVfilt{indd,1}])); % nº RyRs event
end

end