function []=guardaCSVs(Qcom, folder, RF, th, fr, tn, roiR, thI, minR, maxR, sr, as)

ind=find(folder=='\',1,'last');
filename=folder(ind+1:end);

bons={};
cb=0;

mals={};
cm=0;

% Capçalera
% Bons
cb=cb+1;bons{cb,1}=['Folder:;', folder,';'];
cb=cb+1;bons{cb,1}=['Global parameters:;'];
cb=cb+1;bons{cb,1}=['Number of images:;',num2str(tn),';'];
cb=cb+1;bons{cb,1}=['Number of frames:;',fr,';'];
cb=cb+1;bons{cb,1}=['Global threshold:;',num2str(th),';'];
cb=cb+1;bons{cb,1}='';
cb=cb+1;bons{cb,1}=['Filtering parameters:;'];
cb=cb+1;bons{cb,1}=['ROI radius:;',num2str(roiR),';'];
cb=cb+1;bons{cb,1}=['I threshold:;',num2str(thI),';'];
cb=cb+1;bons{cb,1}=['Min radius:;',num2str(minR),';'];
cb=cb+1;bons{cb,1}=['Max radius:;',num2str(maxR),';'];
cb=cb+1;bons{cb,1}='';

% Dolents
mals=bons; cm=cb;

% gruix de ryrs
for ii=1:length(Qcom)
    Qi=Qcom{ii};
    filenamei=[filename '_' num2str(ii)];
    
    % Imatge i paràmetres
    cb=cb+1;bons{cb,1}='';cb=cb+1;bons{cb,1}=['Im' num2str(ii) ';'];
    cb=cb+1;bons{cb,1}=['Zline angle:;',num2str(as(ii)),';'];
    cb=cb+1;bons{cb,1}=[';RyR index; X ; Y; Mean intensity; Max intensity; Radius; Dist NN;Zline; Layer'];
   % if sr==1,
        cm=cm+1;mals{cm,1}='';cm=cm+1;mals{cm,1}=['Im' num2str(ii) ';'];
        cm=cm+1;mals{cm,1}=[';RyR index; X ; Y; Mean intensity; Max intensity; Radius; Dist NN;Zline; Layer'];
    %end
    
    % Mesures de cada RyRs
    for kk=2:size(Qi,1)
        if(strcmp(Qi{kk,9},'Ok')),
            cb=cb+1;bons{cb,1}=[';' num2str(kk-1) ';' num2str(Qi{kk,1}) ';'...
                num2str(Qi{kk,2}) ';' num2str(Qi{kk,3}) ';'...
                num2str(Qi{kk,4}) ';' num2str(Qi{kk,5}) ';'...
                num2str(Qi{kk,6}) ';' num2str(Qi{kk,7}) ';'...
                num2str(Qi{kk,8}) ';'];
        else
            cm=cm+1;mals{cm,1}=[';' num2str(kk-1) ';' num2str(Qi{kk,1}) ';'...
                num2str(Qi{kk,2}) ';' num2str(Qi{kk,3}) ';'...
                num2str(Qi{kk,4}) ';' num2str(Qi{kk,5}) ';'...
                num2str(Qi{kk,6}) ';' num2str(Qi{kk,7}) ';'...
                num2str(Qi{kk,8}) ';'];
        end
    end
end

% muntaCSV(file,separator,variables)
muntaCSV([folder,'\',RF,'\RyRs_', filenamei '.csv'],';',bons);
muntaCSV([folder,'\',RF,'\RyRsRejected_', filenamei, '.csv'],';',mals);

end