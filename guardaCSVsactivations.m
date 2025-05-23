function []=guardaCSVsactivations(EVfilt,folder,RF,ampF,tauF,rorF,t2pF,dhmF,factamp)

ind=find(folder=='\',1,'last');
filename=folder(ind+1:end);

bons={};
cb=0;

mals={}; tot={};

% Capçalera
% Bons

cb=cb+1;bons{cb,1}=['Global parameters:;'];
cb=cb+1;bons{cb,1}='';
cb=cb+1;bons{cb,1}=['Events filtering parameters:;'];
cb=cb+1;bons{cb,1}=['Min amplitude:;',num2str(ampF),';;'];
cb=cb+1;bons{cb,1}=['Min tau:;',num2str(tauF(1)),';Max tau:;',num2str(tauF(2)),';;'];
cb=cb+1;bons{cb,1}=['Min RoR:;',num2str(rorF(1)),';Max RoR:;',num2str(rorF(2)),';;'];
cb=cb+1;bons{cb,1}=['Min t2p:;',num2str(t2pF(1)),';Max t2p:;',num2str(t2pF(2)),';;'];
cb=cb+1;bons{cb,1}=['Min DHM:;',num2str(dhmF(1)),';Max DHM:;',num2str(dhmF(2)),';;'];
cb=cb+1;bons{cb,1}='';

% Dolents
mals=bons; tots=bons; cm=cb; cc=cb;
cb=cb+1;bons{cb,1}=[';ROI; ROI index; X; Y; t; amp; tau; RoR; t2p; DHM; AMP; r2;'];
cm=cm+1;mals{cm,1}=[';ROI;ROI index; X; Y; t; amp; tau; RoR; t2p; DHM; AMP; r2; Note Rejected;'];
cc=cc+1;tot{cc,1}=[';ROI;ROI index; X; Y; t; amp; tau; RoR; t2p; DHM; AMP; r2; Note;'];

% Mesures de cada RyRs
for kk=1:size(EVfilt,1)
    if(strcmp(EVfilt{kk,20},'Ok')),
        cb=cb+1;bons{cb,1}=[';' num2str(EVfilt{kk,1}) ';'...
            num2str(EVfilt{kk,2}) ';' num2str(EVfilt{kk,4}) ';'...
            num2str(EVfilt{kk,5}) ';' num2str(EVfilt{kk,6}) ';'...
            num2str(EVfilt{kk,11}) ';' num2str(EVfilt{kk,16}) ';'...
            num2str(EVfilt{kk,13}) ';' num2str(EVfilt{kk,14}) ';'...
            num2str(EVfilt{kk,15}) ';' num2str(EVfilt{kk,12}) ';'...
            num2str(EVfilt{kk,17}) ';'];
    else
        cm=cm+1;mals{cm,1}=[';' num2str(EVfilt{kk,1}) ';'...
              num2str(EVfilt{kk,2}) ';' num2str(EVfilt{kk,4}) ';'...
            num2str(EVfilt{kk,5}) ';' num2str(EVfilt{kk,6}) ';'...
            num2str(EVfilt{kk,11}) ';' num2str(EVfilt{kk,16}) ';'...
            num2str(EVfilt{kk,13}) ';' num2str(EVfilt{kk,14}) ';'...
            num2str(EVfilt{kk,15}) ';' num2str(EVfilt{kk,12}) ';'...
            num2str(EVfilt{kk,17}) ';' num2str(EVfilt{kk,20}) ';'];
    end
    
            cc=cc+1;tot{cc,1}=[';' num2str(EVfilt{kk,1}) ';'...
             num2str(EVfilt{kk,2}) ';' num2str(EVfilt{kk,4}) ';'...
            num2str(EVfilt{kk,5}) ';' num2str(EVfilt{kk,6}) ';'...
            num2str(EVfilt{kk,11}) ';' num2str(EVfilt{kk,16}) ';'...
            num2str(EVfilt{kk,13}) ';' num2str(EVfilt{kk,14}) ';'...
            num2str(EVfilt{kk,15}) ';' num2str(EVfilt{kk,12}) ';'...
            num2str(EVfilt{kk,17}) ';' num2str(EVfilt{kk,20}) ';'];
end


% muntaCSV(file,separator,variables)
muntaCSV([folder,'\',RF,'\Events_filt_' num2str(factamp) '_',filename,'.csv'],';',bons);
%muntaCSV([folder,'\',RF,'\EventsRejected_',filename,'.csv'],';',mals);
muntaCSV([folder,'\',RF,'\Events_nofilt_' num2str(factamp) '_',filename,'.csv'],';',tot);

end