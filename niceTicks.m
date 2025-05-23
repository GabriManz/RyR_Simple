function [tiks] = niceTicks(yl)
% posicions dels ticks per un graph on els ylims son yl

nombre=[2 10];% interval de nombre de ticks q volem

aux=[1 2 5];% primera xifra sig dels tiks q poden sortir
%complet=zeros(1,90);
L=length(aux);
for ii=1:30% recorrem rang de potencies
k=10^(ii-15);

%complet = [complet aux*k];
complet(((ii-1)*L)+1:((ii-1)*L)+L)=aux*k;
end




interv=yl(2)-yl(1);

nticks=interv./complet;

aux=find((nticks>=nombre(1))&(nticks<=nombre(2)));

intervals=complet(aux(1));% ens qdem el q mes



tiks=floor(yl(1)):intervals:ceil(yl(2));% generem rang mes gran q ylim
%tiks=tiks-mod(tiks,intervals);% anys despres ho trec pq no se q fa
tiks(tiks<yl(1))=[];
tiks(tiks>yl(2))=[];




end