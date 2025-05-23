function [EVcorr]=EventParameters(EV,folder,RF,nnr,signals,plt,DT,DX,wt,pos,V2,events)
% Calcula els paràmetres dels events de la matriu.
% Representa els paràmetres calculats (plt==1)


% Prepara matriu on es guarda la info dels paràmetres:
EVcorr=[EV,zeros(size(EV,1),18)]; EVcorr=num2cell(EVcorr);

% Calcul parametres sparks
extra='prof_spks';
if (exist([folder '\' RF '\' extra],'dir')==0), mkdir([folder '\' RF '\' extra]);end
ensenya('Event parameters. 00%');fprintf('\b');
% plt=0;
% wt=ceil(100/DT); % ms 150
normVal=10;
normVal=max(max(signals));
% cm=[255 255 255;0 0 0;10 193 168;252 168 0;156 0 252]/255; % Colormap
% normVal=[ceil(max(max(signals))*10)/10]; normVal=[ceil(10*normVal/3)/10 normVal]; % Escala --> 5
for ii=1:size(EVcorr,1)
    numeret=round((ii/size(EVcorr,1))*100);
    if(numeret<10),numeret=['0' num2str(numeret)];else numeret=num2str(numeret);end
    fprintf(['\b\b\b\b ' numeret '%%']);
    ts=EVcorr{ii,6}; sig=signals(EVcorr{ii,1},:);
    ind1=max([1,ts-wt]);ind2=min([nnr,ts+3*wt]);
    tv=[ind1:ind2];[~,centre]=intersect(tv,ts);
    [sortida r2 I bline C] = sparkParametersFast2([500,500],tv*DT,sig(tv),centre,plt,['Event ' num2str(ii) '; ROI ' num2str(EVcorr{ii,1})],length(tv),normVal,superjet(6,'wktov5'),DT,nnr);
    % sortida=[amp,tau,ror,t2p,dhm,AMP];
    if(plt==1),imwrite(I,[folder '\' RF '\' extra '\Evt_' num2str(ii) '.png']);end
    EVcorr{ii,6}=tv(C); % t
    EVcorr{ii,11}=sortida(1); % amp
    EVcorr{ii,12}=sortida(6); % AMP
    EVcorr{ii,13}=sortida(3); % ror
    EVcorr{ii,14}=sortida(4); % t2p
    EVcorr{ii,15}=sortida(5); % dhm
    EVcorr{ii,16}=sortida(2); % tau
    EVcorr{ii,17}=r2; EVcorr{ii,18}=bline;
    EVcorr{ii,21}=sortida(7);  EVcorr{ii,22}=sortida(8);  EVcorr{ii,23}=sortida(9);  EVcorr{ii,24}=sortida(10); % FD 25,50,75,90
end
fprintf(' ');
save([folder,'\',RF,'\detectionParam.mat'],'EVcorr','wt','DT','DX','events','signals','V2','pos','EV');

end