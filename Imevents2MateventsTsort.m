function [EV]=Imevents2MateventsTsort(events,pos,rr,Q)
% Passa de la imatge de la detecció d'events en les senyals a una matriu
% amb tots els events:
%   * En la matriu hi ha info sobre de quina ROI provenen, RyR pos i temps
% Finalment els ordena en funció del temps.

% Matriu events
% ROI ROIryr x y t parametresSpk cluster
% RyRs ROI 1:rrr, ROI>rrr no RyRs
ind=find(events>0);
if iscell(Q), Q=Q{1}; end;
numRyRm=max(Q(:,1)); % Mira el num maxim assignat a un RyR
EV=zeros([length(ind),6]);
rrr=rr;
ct=1; ct2=1;
for ii=1:rr
    inde=find(events(ii,:)==1); % Passa per a cada ROI mirant si té events la senyal
    for jj=1:length(inde)
        if ii<=rrr, ryr=Q(ii,1); EV(ct,3)=1; % RyR 1/no RyR 0
        else ryr=numRyRm+ct2; ct2=ct2+1;end % Si és RyR o no RyR assignar valor diferent
        EV(ct,1)=ii; % Num ROI senyal
        EV(ct,2)=ryr; % Num ROI ryr
        EV(ct,4)=pos(ii,2); % x
        EV(ct,5)=pos(ii,1); % y
        EV(ct,6)=inde(jj); % Frame
        ct=ct+1;
    end
end

% tt=zeros(size(EV(:,6)));
[~,I] = sort(EV(:,6));
EV2=zeros(size(EV));
for ii=1:length(I), EV2(ii,:)=EV(I(ii),:);end
EV=EV2;
clear EV2
        
end