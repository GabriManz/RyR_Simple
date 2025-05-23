function [EVcorr,blS,blSstd]=FilterEvents(folder,RF,ampF,tauF,rorF,t2pF,dhmF,r2F,EVcorr,ev,signals,eventsW,zzRyR,Iryr,pos,DX,DT,rr,factamp,nnr,limmad,lim)
% Filtra els events segons els paràmetres estimats.

cc=1;
ampFo=ampF;
[blS,blSstd]=BlandBlstd(ev,eventsW,DT,rr,signals,nnr); % Bl i Sigma bl
limmad2=(limmad./lim)*(1/3); 
% limmad2=(limmad./lim);% Normalització de MAD(Fo) --> MAD(Fo)norm=MAD(Fo)/Fo 
% --> Entre això i el factamp acaba quedant --> ampF = MAD(Fo)norm * 2.5
% --> Amp ha de ser 2.5 vegades més gran que el soroll!
[pos2, pos, ind]=EliminaRyRsFronteres(zzRyR,Iryr,pos,DX); % Out RyRs dels bordes!
factampo=6;
for ii=1:size(EVcorr,1)
    ampF=ampFo;
    ampM=(factamp/3)*blSstd(EVcorr{ii,1});%+blS(EVcorr{ii,1});
    %if ampF<ampM,ampF=ampM; end
    ampMM=limmad2(EVcorr{ii,1})*factamp*1.2;% factamp=6;
     % ampMM=limmad2(EVcorr{ii,1})*factamp;%*1.2; factamp=2.4
    ampF=ampMM;
    if ~isnan(EVcorr{ii,11}), EVcorr{ii,20}='Ok';
        if EVcorr{ii,12}>ampF,  EVcorr{ii,20}='Ok'; % AMP (12),  no amp (11) (AMP-bl)
            if EVcorr{ii,16}>tauF(1),   EVcorr{ii,20}='Ok';
                if EVcorr{ii,16}<tauF(2),   EVcorr{ii,20}='Ok';
                    if EVcorr{ii,13}>rorF(1),   EVcorr{ii,20}='Ok';
                        if EVcorr{ii,13}<rorF(2),   EVcorr{ii,20}='Ok';
                            if EVcorr{ii,14}>t2pF(1),   EVcorr{ii,20}='Ok';
                                if EVcorr{ii,14}<t2pF(2),   EVcorr{ii,20}='Ok';
                                    if EVcorr{ii,15}>dhmF(1),   EVcorr{ii,20}='Ok';
                                        if EVcorr{ii,15}<dhmF(2),  EVcorr{ii,20}='Ok';
                                            if EVcorr{ii,17}>r2F,  EVcorr{ii,20}='Ok';
                                                if all(EVcorr{ii,1}~=ind), EVcorr{ii,20}='Ok';
                                                    EVcorr{ii,19}=1;
                                                    % EVfilt{cc,:}=EVcorr2{ii,:}; cc=cc+1; % Supera els filtres!
                                                else
                                                    EVcorr{ii,19}=0; EVcorr{ii,20}='RyR in the border';
                                                end
                                            else
                                                EVcorr{ii,19}=0; EVcorr{ii,20}='r2 too low';
                                            end
                                        else
                                            EVcorr{ii,19}=0; EVcorr{ii,20}='DHM too big';
                                        end
                                    else
                                        EVcorr{ii,19}=0; EVcorr{ii,20}='DHM too small';
                                    end
                                else
                                    EVcorr{ii,19}=0; EVcorr{ii,20}='t2p too big';
                                end
                            else
                                EVcorr{ii,19}=0; EVcorr{ii,20}='t2p too small';
                            end
                        else
                            EVcorr{ii,19}=0; EVcorr{ii,20}='RoR too big';
                        end
                    else
                        EVcorr{ii,19}=0; EVcorr{ii,20}='RoR too small';
                    end
                else
                    EVcorr{ii,19}=0; EVcorr{ii,20}='tau too big';
                end
            else
                EVcorr{ii,19}=0; EVcorr{ii,20}='tau too small';
            end
        else
            EVcorr{ii,19}=0; EVcorr{ii,20}=['Low amp ' num2str(ampF)];
        end
    else
        EVcorr{ii,19}=0; EVcorr{ii,20}='Is NaN';
    end
    ampF=ampFo;
end

%         [a,b]=size(Iryr);
%         IQ=zeros([a,b]);
%         b=[EVcorr{:,4}];a=[EVcorr{:,5}];
%         for ii=1:length(a)
%             IQ(a(ii),b(ii))=1;
%         end
%         [~,I] = sort([EVcorr{:,6}]);  EV2=zeros(size(EV));for ii=1:length(I), for jj=1:size(EVcorr,2), EV2{ii,jj}=EVcorr{I(ii),jj};end;end
guardaCSVsactivations(EVcorr,folder,RF,ampF,tauF,rorF,t2pF,dhmF,factampo);

end