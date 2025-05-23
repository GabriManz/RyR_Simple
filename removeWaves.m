function [EVfilt2,EVfilt]=removeWaves(EVfilt,EVfilt2)

indw=find([EVfilt{:,25}]==1); % Busca waves
if ~isempty(indw),
    inds=[1:size(EVfilt,1)]; inds(intersect(inds,indw))=[];
    % Busca ets en waves i en els sparks
    etw=[EVfilt{indw,10}];
    ets=[EVfilt{inds,10}];
    indd=intersect(ets,etw); % Mira et que coincideixen
    if ~isempty(indd),
        inde=[];
        for jj=1:length(indd), indei=find(ets==indd(jj)); inde=[inde , indei]; end
        for jj=1:length(inde),inds(inds==inde(jj))=[];end % Elimina sparks amb la mateix et que waves
    end
    EVfilt2=EVfilt;
    EVfilt=cell([length(inds), size(EVfilt,2)]);
    for jj=1:length(inds), for kk=1:size(EVfilt2,2), EVfilt{jj,kk}=EVfilt2{inds(jj),kk}; end, end
end

end