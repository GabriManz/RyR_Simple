function [score] = comparaTextos(s1,s2)


l1=length(s1);
l2=length(s2);

if(l1>l2),aux=s1;s1=s2;s2=aux;l1=length(s1);l2=length(s2);end% volem s1 la curta

score=0;
for ii=l1:-1:1 % per cada possible longitud de la curta (desde la seq completa fins a caracter individual)
sc=0;
    for jj=1:l1-ii+1 % per cada possible encaix de la longitud ii en la curta
        scc=0;
        scc1=0;
      fin=jj-1;
   ori=jj+ii;
   kk=findstr(s1(jj:jj+ii-1),s2);%busca la sequencia amb longitud ii de la curta dins la string llarga

   if(~isempty(kk))% si ha trobat la sequencia exacte
   scc=scc+ii;% sumem un score acumulat igual a la longitud de sequencia encaixada
   scc1=max([scc1 ii]);% el score maxim per aquest encaix 
   end
   
   kk=findstr(s1(1:jj-1),s2(1:fin)); %busca el tram anterior a lencaix de la curta dins la string llarga
   if(~isempty(kk))
   scc=scc+jj-1;% sumem un score acumulat igual a la longitud del tram anterior a la sequencia encaixada
   scc1=max([scc1 jj-1]);% el score maxim per aquest encaix 
   end    
   
   kk=findstr(s1(jj+ii:end),s2(ori:end)); %busca el tram posterior a lencaix de la curta dins la string llarga
   if(~isempty(kk))
   scc=scc+l1-jj-ii+1;% sumem un score acumulat igual a la longitud del tram posterior a la sequencia encaixada
   scc1=max([scc1 l1-jj-ii+1]);% el score maxim per aquest encaix 
   end
   sc=max([sc (scc+scc1)*.5]);% SCORE DE LENCAIX ES LA MITJA ENTRE EL MAXIM I LACUMULAT
    end
   score=max([score sc]);% el SCORE absolut es el maxim dels obtinguts en totes les longituds
   
end


score=score/l2;% puntuacio normalitzada a la llarga

end