function [allx,ally] = puntsEnmig(xx,yy)
%troba els punts entre la parella de punts donada


if((xx(1)==xx(2))&&(yy(1)==yy(2)))
allx=xx;ally=yy;
else
allx=[];ally=[];



rangy=yy(1):((yy(2)-yy(1))/abs((yy(2)-yy(1)))):yy(2);
rangx=xx(1):((xx(2)-xx(1))/abs((xx(2)-xx(1)))):xx(2);

N=max(length(rangx),length(rangy));

rangy=(linspace(yy(1),yy(2),N));
rangx=(linspace(xx(1),xx(2),N));

[m,n] = punts2recta(xx(1),yy(1),xx(2),yy(2));
if(abs(m)==Inf),
    [m,n] = punts2recta(yy(1),xx(1),yy(2),xx(2));
    for kk=rangy
        allx=[allx (m*kk+n)];
        ally=[ally kk];
    end
else
    if(m==0)
        for kk=rangx
            Y=round(m*kk+n);
            ally=[ally (m*kk+n)];
            allx=[allx kk];
        end
    else
        for kk=rangy
            allx=[allx ((kk-n)/m)];
            ally=[ally kk];
        end
    end
end


allx=round(allx);
ally=round(ally);

        
end  

    
    
end