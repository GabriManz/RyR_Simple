function [maxims,minims] = findPeaks6(data_in)
%Detect maxia

        Deriv=diff(data_in);

        P=Deriv(1:end-1).*Deriv(2:end);

        candidats=find(P<0)+1;

        maxims=candidats(Deriv(candidats)<0);

        minims=candidats(Deriv(candidats)>0);
end
 
