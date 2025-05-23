function [h]=subNM(n,m,q,marg)
% subplot fix margin, creates axes in current figure such that it is the
% q-th out of a n-by-m grid with no margins in between. 
if(nargin<4),marg=0;end

fila=ceil(q/m);
columna=q-((fila-1)*m);


pos(1)=min((columna*marg)+(columna-1)*(1-(m+1)*marg)/m);
pos(2)=min(-((1-(n+1)*marg)/n)+1-((fila*marg)+(fila-1)*(1-(n+1)*marg)/n));
pos(3)=(columna(end)-columna(1)+1)*((1-(m+1)*marg)/m)+(columna(end)-columna(1))*marg;
pos(4)=(fila(end)-fila(1)+1)*((1-(n+1)*marg)/n)+(fila(end)-fila(1))*marg;


h=axes('units','normalized','position',pos);
    




end