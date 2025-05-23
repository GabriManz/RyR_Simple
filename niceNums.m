function [M] = niceNums(N,opt)
if(nargin==1),opt=2;end
M=num2str(N);

if(~isnan(N))



e=strfind(lower(M),'e');

if(~isempty(e))
    
  xp=M(e+1:end);
  M=M(1:e-1);
  
   %opt=0; 
end

MM=M;
p=strfind(M,'.');
if(~isempty(p))
    if(opt==0)
    M=num2str(round(str2num(M)));
    end
    
    if(opt>0)
        M=M(1:p);
        extr='';
        for kk=1:length(MM(p+1:end))
           if(MM(p+kk)=='0'),extr=[extr '0'];else break;end 
        end
        decimals=num2str(round((10^opt)*str2num(['0.' MM(p+1:end)])));
        decimals=[extr decimals];
        if(length(decimals)>opt)
        M=num2str(str2num(M)+(str2num(decimals)/(10^opt)));
        else
        M=[M decimals];
        end
    end

  
end




if(~isempty(e))
   M=[M 'e' xp]; 
end

end

end