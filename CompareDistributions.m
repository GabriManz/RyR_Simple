function [] = CompareDistributions(varargin)
%
%  CompareDistributions(var); where var is cell of distribution values.
%
%  CompareDistribution(var,'optionName',optionvalue,...);
%  Allows specific options amongst the following:
%
% 'tst' test to apply either 'kstest', 'ranksum', or 'ttest'(default)
% 'dsp' display type either 'bar' or 'box',
% 'lab' labels for each set in var,
% 'ori' orientation of plot either 'vert' or 'hori',
% 'tit' plot title,
% 'uni' units or text for the axis label,
% 'sig' significance shown as either '*' or 'p' to display p-values,
% 'fig' call new figure either 1 or 0,
% 'prt' print file either '0' or a path to the file,
% 'fs'  fontsize in figure,
% 'dc'  datacolors can either be a single color (1 by 3) or a matrix of colors (N by 3),
% 'ec'  errorcolors can either be a single color (1 by 3) or a matrix of colors (N by 3),

argind=0;
argind=argind+1;if(nargin<argind),help CompareDistributions;error('No samples!');end
var=varargin{1};
if(mod(length(varargin),2)==0),help CompareDistributions;error('Options must be called by couples of parameters; propertiy name as string followed by variable');end
tst='ttest';
dsp='box';
lab=cell(length(var));for ii=1:length(lab),lab{ii}=num2str(ii);end;
ori='hori';
FS=11;
pv='*';
titol='';
varuni='';
prt='0';
avisfinal=' ';
nf=1;
datacolors=[0 0 .6];
errorcolors=[.6 0 0];
wid=0.3;
ii=0;
nbf=1; % Bonferroni correction
for jj=2:2:length(varargin)
    % ii=ii+2;
    ii=jj;
    if(ischar(varargin{ii}))
        switch lower(varargin{ii})
            case 'tst'
                tst=varargin{ii+1};
            case 'dsp'
                dsp=varargin{ii+1};
            case 'lab'
                lab=varargin{ii+1};
            case {'ori'}
                ori=varargin{ii+1};
            case {'fs','font'}
                FS=varargin{ii+1};
            case 'tit'
                titol=varargin{ii+1};
            case 'uni'
                varuni=varargin{ii+1};
            case 'fig'
                nf=varargin{ii+1};
            case 'sig'
                pv=varargin{ii+1};
            case 'prt'
                prt=varargin{ii+1};
            case {'dc'}
                datacolors=varargin{ii+1};
            case 'ec'
                errorcolors=varargin{ii+1};
            case 'bf'
                nbf=varargin{ii+1};
            case 'wi'
                wid=varargin{ii+1};
            otherwise
                avisfinal=varargin{ii};
        end
    else
        ii=ii+1;
    end
    
end



labs=lab;
N=length(labs);
if(size(datacolors,1)<N),for ii=1:N,datacolors(ii,:)=datacolors(1,:);end;end
if(size(errorcolors,1)<N),for ii=1:N,errorcolors(ii,:)=errorcolors(1,:);end;end
if(strcmp(tst,'ttest')==1),ts=1;end
if(strcmp(tst,'kstest')==1),ts=2;end
if(strcmp(tst,'ranksum')==1),ts=3;end
if(strcmp(tst,'kwtest')==1),ts=4;end
if(ori=='hori'),oo=1;end
if(ori=='vert'),oo=2;end
if(strcmp(prt,'0')),imp=0;else imp=1;end
%

testos={'t-test','Kolmogorov-Smirnov test','Wilcoxon rank sum test','KRUSKALWALLIS'};
orient={'horizontal','vertical'};
ch='rgbcymkw';
co=[1,0,0;0 1 0;0 0 1;0 1 1;1 1 0;1 0 1;0 0 0;1 1 1];
if(dsp=='box')
    for ii=1:N,
        for jj=1:size(co,1),
            d(jj)=norm(co(jj,:)-errorcolors(1,:));
        end
        [aux,quin]=min(d);
        dst(ii)=aux;
        qui(ii)=quin;
    end
    [aux,quin]=min(dst);
    errcol=ch(qui(quin));
end
unitats=varuni;

V=[];C=[];
for ii=1:N
    v=var{ii};v(isnan(v))=[];
    if(size(v,1)>size(v,2)),v=v';end
    var{ii}=v;
    mitj(ii)=mean(v);
    medi(ii)=median(v);
    %stds(ii)=std(v);
    % CARMETA:
    stds(ii)=std(v)/sqrt(length(v));
    %
    V=[V,v];
    C=[C,ii*ones(1,length(v))];
    q9(ii)=quantile(v,.95);
end

if(nf),figure;else gca;end
hold on;
desp=.05;
correccio=0;
if(dsp=='box')
    %for ii=1:N
    
    boxplot(V,C,'orientation',orient{oo},'colors',datacolors,'symbol',['+' errcol]);
    %boxplot(var{ii},'positions',ii,'colors',datacolors,'symbol',errcol(ii));
    %hold on;
    %    end
    %if(oo==1),correccio=.15;else correccio=-.45;end
    if(oo==1),correccio=.15;else correccio=-.25;end
    % correccio=-.2
end
for ii=1:N
    if(oo==1)
        h=text(medi(ii),ii-.2-correccio,[' (' num2str(length(var{ii})) ')']);
    else
        h=text(ii-correccio,mitj(ii),[' (' num2str(length(var{ii})) ')']);
        set(h,'verticalalignment','bot');
    end
    if(dsp=='bar')
        if(oo==1)
            set(h,'horizontalalignment','left','verticalalignment','middle');
            h=line([mitj(ii) mitj(ii)+stds(ii)],[ii ii]);set(h,'color',errorcolors(ii,:));
            h=line([mitj(ii)+stds(ii) mitj(ii)+stds(ii)],[ii-desp ii+desp]);set(h,'color',errorcolors(ii,:));
        else
            set(h,'horizontalalignment','left','verticalalignment','bottom');
            h=line([ii ii],[mitj(ii) mitj(ii)+stds(ii)]);set(h,'color',errorcolors(ii,:));
            h=line([ii-desp ii+desp],[mitj(ii)+stds(ii) mitj(ii)+stds(ii)]);set(h,'color',errorcolors(ii,:));
        end
    else
        set(h,'horizontalalignment','center','verticalalignment','middle');
    end
end
hold on;
if(dsp=='bar')
    for ii=1:N
        if(oo==1)
            barh(ii,mitj(ii),'facecolor',datacolors(ii,:),'edgecolor',datacolors(ii,:));
        else
            bar(ii,mitj(ii),wid,'facecolor',datacolors(ii,:),'edgecolor',datacolors(ii,:));
        end
    end
end

% UNITATS
if(oo==1)
    set(gca,'ytick',[1:N],'yticklabel',labs);
    xlabel(unitats,'fontsize',FS-1);
else
    set(gca,'xtick',[1:N],'xticklabel',labs);
    ylabel(unitats,'fontsize',FS-1);
end
if(dsp=='bar'),oc0=max(mitj+stds);end
if(dsp=='box'),oc0=max(q9);end
if(imp==1),pd=-.03;fd=.04;ocF=1.5;ma=.2;else pd=-.06;fd=.05;ocF=2;ma=.5;end
fd=.1; %%%%%%%%%%%%%%%%%%%%%% MOVIMENT ESTRELLETES!!!
d2=oc0*fd;d3=.05;
oc0=oc0+d2;or=oc0;
ocu=zeros(N,max(0,nchoosek(N,2)+1));

disp(' ');
% if(strcmp(titol,''))
% titol=['Analysis using ' testos{ts} '.'];
% else
% titol=['Analysis for ' titol ' using ' testos{ts} '.'];
% end
disp(titol);
title(titol,'fontsize',FS+2);
szc=10;

% BONFERRONI
pbf1=.05/nbf;
pbf2=.01/nbf;

% TEST
for ii=1:N-1
    FA=labs{ii};
    for jj=ii+1:N
        FAA=labs{jj};
        try,
            if(ts==1),[H,p]=ttest2(var{ii},var{jj});end
            if(ts==2),[H,p]=kstest2(var{ii},var{jj});end
            if(ts==3),p=ranksum(var{ii},var{jj});end
            %if(ts==4),p=kruskalwallis(var{ii},var{jj});end
            disp([FA ' vs ' FAA]);
            disp(['p-value: ' num2str(p)]);
            if(p<pbf1)
                if(pv=='*')
                    sen=' *';
                    % if(p<pbf1),sen=[sen '*'];end
                    if(p<.01),sen=[sen '*'];end
                    if(p<.001),sen=[sen '*'];1
                    end
                     % if(p<.0001),sen=[sen '***'];end
                     %if(p<.0001),sen=[sen '****'];end
                else
                    sen=[' p=' niceNums(p)];
                    szc=-1;
                end
                espai=0;or=oc0;
                while(espai==0)
                    if(sum(ocu(ii:jj,1+round((or-oc0)/(ocF*d2))))<=1)
                        espai=1;
                    else
                        or=or+ocF*d2;
                    end
                end
                if(oo==1)
                    h=line([or or+d2],[ii+d3 ii+d3]);set(h,'color','k');
                    h=line([or or+d2],[jj-d3 jj-d3]);set(h,'color','k');
                    h=line([or+d2 or+d2],[ii+d3 jj-d3]);set(h,'color','k');
                    h=text(or+d2,pd+(ii+jj)*.5,sen);set(h,'fontsize',FS+szc,'verticalalignment','middle','horizontalalignment','left');
                else
                    h=line([ii+d3 ii+d3],[or or+d2]);set(h,'color','k');
                    h=line([jj-d3 jj-d3],[or or+d2]);set(h,'color','k');
                    h=line([ii+d3 jj-d3],[or+d2 or+d2]);set(h,'color','k');
                    h=text((ii+jj)*.5,or+d2*.2,sen);set(h,'fontsize',FS+szc,'horizontalalignment','center','verticalalignment','bottom');
                end
                ocu(ii:jj,1+round((or-oc0)/(ocF*d2)))=1;
                
            end
        end
    end
end

if(oo==1)
    kk=xlim();
    xlim([kk(1) oc0+ocF*d2*(find(sum(ocu)==0,1,'first')-ma)]);
else
    kk=ylim();
    ylim([kk(1) oc0+ocF*d2*(find(sum(ocu)==0,1,'first')-ma)]);
end
set(gca,'fontsize',FS);
%ocu=ocu(:,1:find(sum(ocu)==0,1,'first'));
%figure;imagesc(ocu(end:-1:1,:));
if(strcmp(avisfinal,' ')==0)
    
    disp(' ');
    disp('warning');
    disp(['uncomprehesible parameter: ''' avisfinal '''.']);
    
end
disp(' ');

end