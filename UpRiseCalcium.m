function [EVfilt]=UpRiseCalcium(EVfilt, DT,signals)

            wps=round(25/DT);
            wbl=round(150/DT);
            wpost=round(300/DT);
            for ii=1:size(EVfilt,1)
                t=[EVfilt{ii,6}];
                sig=signals(EVfilt{ii,1},:);
                
                in=max([1 t-wbl]);v=[in:t];
                sw=sig(v);
                if length(in:t)>6,
                    quarticMA = sgolayfilt(sig, 8, 19);
                    swf=quarticMA(in:t); % https://es.wikipedia.org/wiki/Filtro_de_Savitzky%E2%80%93Golay
                    pend=zeros(1,length(sw)-1); pendf=pend;
                    for jj=2:length(sw)-1
                        pend(jj)=(sw(jj)-sw(1))/jj;
                        pendf(jj)=(swf(jj)-swf(1))/jj;
                    end
                    [MM,mm]=findPeaks6(pendf);
                    if ~isempty(mm),
                        inblsw=mm(end); inbl=v(inblsw)+1;
                    else
                        [~,inbl]=min(sw);
                    end
                    
                else
                    [~,inbl]=min(sw);
                end
                EVfilt{ii,23}=inbl;
            end
            
end