function [Q2,Qgood,Qr] = filtraryrs(Q, Inorm, thi, minR, maxR)


% Q=[X,Y,IntMitja,IntMax,Rx,Ry]
Qi2=[];
Qi3=[];Qi4=[];

for ii=1:size(Inorm,3)
    
    Qi=Q{ii};
    Qi2=[Qi zeros(size(Qi,1),1)]; % Columna amb info del filtratge
    Qi2=[zeros(1,size(Qi2,2));Qi2];
    Qi2=num2cell(Qi2); % converteix en cell per introduir els comentaris
    
    % Afegim titol
    Qi2{1,1}='X';
    Qi2{1,2}='Y';
    Qi2{1,3}='Mean intensity';
    Qi2{1,4}='Max intensity';
    Qi2{1,5}= 'Radius';
    Qi2{1,6}='Comments';
    
    Qi3=zeros(size(Qi,1),size(Qi,2)+1);
    Qi4=zeros(size(Qi,1),size(Qi,2)+1);
    
    cont=1; cont2=1;
    for kk=1:length(Qi)
        if (Qi(kk,4)>thi),
            if (Qi(kk,5)>minR),
                if (Qi(kk,5)<maxR),
                    Qi2{kk+1,9}='Ok';
                    Qi3(cont,:)=[kk Qi(kk,:)]; cont=cont+1;
                else
                    Qi2{kk+1,9}='Radius too large';
                    Qi4(cont2,:)=[kk Qi(kk,:)]; cont2=cont2+1;
                end
            else
                Qi2{kk+1,9}='Radius too small';
                Qi4(cont2,:)=[kk Qi(kk,:)]; cont2=cont2+1;
            end
        else
            Qi2{kk+1,9}='Low intensity';
             Qi4(cont2,:)=[kk Qi(kk,:)]; cont2=cont2+1;
        end
    end
    Q2{ii}=Qi2;
    Qgood{ii}=Qi3(1:cont-1,:);
    Qr{ii}=Qi4(1:cont2-1,:);
end




end

