function [Q2, IZ, pealing] = mesures(Isum, Q, radiR, DX, as, IQsum,pg)

[a,b]=size(Isum(:,:,1));
Q2=[];
IZ=zeros(size(Isum));
Iz=zeros(size(Isum(:,:,1)));

for ii=1:size(Isum,3)
    Im=Isum(:,:,ii);IQ=IQsum(:,:,ii);
    angle=as(ii);
    Qi=Q{ii};
    Qi2=Qi;
    [zlines Iz]=clusterRyrs(Qi,DX,Im);
    
    Qi=[Qi zeros(size(Qi,1),6)];
    %%%% Afegir dist to nearest neig
    D=squareform(pdist([Qi(:,1),Qi(:,2)],'euclidean'));D(D==0)=1e9;
    
    for kk=1:size(Qi,1)
        x=Qi(kk,1); y=Qi(kk,2);
        minx=max([1 x-radiR]);
        maxx=min([x+radiR b]);
        miny=max([1 y-radiR]);
        maxy=min([y+radiR a]);
        A=Im(miny:maxy,minx:maxx);
        
        % Mesura intensitat
        Qi(kk,3)=mean(mean(A)); % int mean
        Qi(kk,4)=max(max(A)); % int max
        
        % Mesura radis
        % figure(2);imagesc(A);
        [a2,b2]=size(A);
        
        if (abs(angle)>45)&&(abs(angle)<=135), % Les Zlines estan en l'eix y
            % direccio x
            v1=(1:b2)';
            v2=A(round(b2/2),:)';
            % sigma*2 -> 95'44% de la informació de la distribució
            % sigma -> 68.26% de la informació
            % figure(1);plot(v1,v2);hold on;plot(fx);
        else
            % direccio y
            v1=(1:a2)';
            v2=A(:,round(a2/2));% v2=v2';
            % figure(1);plot(v1,v2);hold on;plot(fy);
        end
        try
            f = fit(v1,v2,'gauss1');
            sigma=f.c1/sqrt(2);
            Qi(kk,5)=sigma*DX; % en um
            %         catch
            %             disp([ii,kk])
        end
        % Distancia
        Qi(kk,6)=min(D(kk,:))*DX; % (um)
        
    end
    
    % Zline de cada RyR
    Qi(:,7)=zlines;
    
    %%%% Fer pealing
    [pealingi, Qp]=pealingRyRs(IQ,DX,pg,Qi2);
    
    Qi(:,8)=Qp(:,3);
    
    Q2{ii}=Qi;
    if size(Isum,3)==1, IZ=Iz; 
        pealing=pealingi; 
    else
        IZ(:,:,ii)=Iz; pealing(:,:,ii)=pealingi;
    end
    
    
end

end

