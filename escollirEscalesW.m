function [inf,sup]=escollirEscalesW2(DT,t,T)
% Escollir escales (small events and big events)

% DT: resolucio (ms/mostres)
% t: temps màxim events petits (ms)
% T: temps màxim events grans (ms)

% Senyal 1 peak:
s=zeros(1,5000);
s(1,2500)=1;

% Transformada:
fin=1000;
w=cwt(s,[1:fin],'gaus2');

gausfamily=zeros(fin,2);
for ii=1:fin
    wi=w(ii,:);
    if ii>100,
        [up,lo] = envelope(wi,round(ii/2),'rms');
        [p,inp]=findpeaks(up,'MinPeakProminence',0.01);
        if any(abs(inp(1:end-1)-inp(2:end))<50),inp=[inp(1),inp(2),inp(end)];end
        [~,c1]=min(up(inp(1):inp(2)));
        [~,c2]=min(up(inp(2):inp(3)));
        ind=inp(1)+c1:inp(2)+c2;
        
    else
        ind=find(wi>0); % peak superior a 0 %canviat 180417
    end
    suporti=length(ind)*DT;
    gausfamily(ii,1)=ii;
    gausfamily(ii,2)=suporti;
end

indinf=find(gausfamily(:,2)>=t,1,'first'); inf=max([1,indinf-1]);
indsup=find(gausfamily(:,2)>=T,1,'first'); sup=max([1,indsup-1]);

if indsup<indinf, indsup=indinf; end;

% plot(s);hold on;plot(w(1,:));

end