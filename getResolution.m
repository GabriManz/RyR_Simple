function [res] = getResolution(ruta)
% pensat per stacks

fid=fopen(ruta,'rt');

S=fscanf(fid,'%s');

fclose(fid);


%%% unitats de longitud en SI
SI{1,1}='Gm';SI{1,2}=1e9;
SI{2,1}='Mm';SI{2,2}=1e6;
SI{3,1}='km';SI{3,2}=1e3;
SI{4,1}='hm';SI{4,2}=1e2;
SI{5,1}='m';SI{5,2}=1e0;
SI{6,1}='dm';SI{6,2}=1e-1;
SI{7,1}='cm';SI{7,2}=1e-2;
SI{8,1}='mm';SI{8,2}=1e-3;
SI{9,1}='µm';SI{9,2}=1e-6;
SI{10,1}='nm';SI{10,2}=1e-9;


res=[0 0];



for kk=1:2
    try
    ii=(kk*2)-1;
    NE=strfind(S,'NumberOfElements');
LE=strfind(S,'Length');
UN=strfind(S,'Unit');
if(((~isempty(NE))&&(~isempty(LE)))&&(~isempty(UN)))

UN=UN(UN>NE(1));

NE=NE(ii);
LE=LE(ii);
UN=UN(ii);

NEC=find(S(NE:end)=='"',2,'first');
LEC=find(S(LE:end)=='"',2,'first');
UNC=find(S(UN:end)=='"',2,'first');

nume=str2double(S(NE+NEC(1):NE+NEC(2)-2));
leng=(S(LE+LEC(1):LE+LEC(2)-2));leng(leng==',')='.';leng=str2num(leng);
leng=abs(leng);
unit=S(UN+UNC(1):UN+UNC(2)-2);



R=0;IN=0;
% busca de quina es tracta
for ii=1:10
    r=comparaTextos(unit,SI{ii,1});
    %disp([num2str(ii) ': ' num2str(r)]);
    if(r>R),R=r;IN=ii;end
end

leng=leng*SI{IN,2}/1e-6;% passem a micres

res(kk)=leng/nume; % per unitat de pixel
%%%




end
    end
end
end