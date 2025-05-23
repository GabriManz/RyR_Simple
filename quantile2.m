function [val] = quantile2(v,f)


v=sort(v);

val= zeros(1,length(f));

for ii=1:length(f)
    ind=round(f(ii)*length(v));
    if(ind<1),ind=1;end
val(ii)=v(ind);
end



end


