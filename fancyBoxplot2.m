function [] = fancyBoxplot2(varargin)

hold on;

argind=0;
argind=argind+1;if(nargin<argind),help CompareDistributions;error('No samples!');end
D=varargin{1};
X=varargin{2};
C=varargin{3};
tst='ttest';
pv='*';
prt='0';
nbf=1;
FS=12;
ii=0;
oo=2;
for jj=3:2:length(varargin)
    ii=ii+2;
    if(ischar(varargin{ii}))
        switch lower(varargin{ii})
            case 'tst'
                tst=varargin{ii+1};
            case 'sig'
                pv=varargin{ii+1};
            case 'prt'
                prt=varargin{ii+1};
                case 'lab'
                lab=varargin{ii+1};
            case 'bf'
                nbf=varargin{ii+1};
        end
    end
end
if(strcmp(prt,'0')),imp=0;else imp=1;end

if(nargin<2),N=1;end
N=unique(X);
if(nargin<3),C=zeros(length(N),3)+.4;end
dst=min(diff(N));
if(nargin<3),C=zeros(length(N),3)+.4;end
dst=min(diff(N));

ms=25;% tamany markers
lw=3;% gruix ralles
sz=.7*dst;% llargada ralla gran
fa=.68;% factor entre ralla mes gran i mes peke
fcL=.4;% factor de color ralles (<1 per enfosquir)
res=50; % nombre de bins per determinar solapament
amp=.3; % factor d'ample maxim dels punts solapats respecte llargada ralla gran
qq=[.1 .5 .9];% quantils a mostrar (hauria de ser simetric)

ed=linspace(min(D),max(D),res);
labs=unique(X);
outl=zeros(1,length(labs));
for jj=1:length(N)
    
    
    
    x=N(jj);
    dist=D(X==x);
    
    
    CP=C(jj,:);
    CL=C(jj,:)*fcL;CL(CL>1)=1;
    
    hi=histc(dist,ed);hh=zeros(size(hi));
    sep=sz*amp/max(hi);
    for ii=1:length(dist)
        on=find(ed<=dist(ii),1,'last');hh(on)=hh(on)+1;
        %hh(on)/hi(on)
        pos=x+(hh(on)*sep-sep*hi(on)*.5)-sep*.5;
        h=plot(pos,dist(ii),'.');   set(h,'markersize',ms,'color',CP);
        h=plot(pos,dist(ii),'ko');set(h,'markersize',ms*.3);% borde negre
    end
    qs=quantile(dist,qq);
    aux=sz*(fa+(1-fa)*linspace(0,1,ceil(length(qs)/2)));aux=[aux(1:end-1) aux(end:-1:1)];
    for ii=1:length(qs)
        h=line([x-aux(ii)/2 x+aux(ii)/2],[qs(ii) qs(ii)]);
        set(h,'linewidth',lw,'color',CL);
    end
    outl(jj)=max(dist);
end
xlim([min(X)-dst/2 max(X)+dst/2]);

% TEST
testos={'t-test','Kolmogorov-Smirnov test','Wilcoxon rank sum test','KRUSKALWALLIS'};
if(strcmp(tst,'ttest')==1),ts=1;end
if(strcmp(tst,'kstest')==1),ts=2;end
if(strcmp(tst,'ranksum')==1),ts=3;end
if(strcmp(tst,'kwtest')==1),ts=4;end

N2=length(labs);
oc0=max(outl);
if(imp==1),pd=-.03;fd=.04;ocF=1.5;ma=.2;else pd=-.06;fd=.05;ocF=2;ma=.5;end
d2=oc0*fd;d3=.05; 
oc0=oc0+d2;or=oc0;
ocu=zeros(N2,max(0,nchoosek(N2,2)+1));
szc=10;

var=cell(1,N2);
for ii=1:N2, var{ii}=[D(X==labs(ii))]; end

pbf1=.05/nbf; % Correcció bonferroni
for ii=1:N2-1
    FA=lab{ii};
    for jj=ii+1:N2
        FAA=lab{jj};
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
                    if(p<pbf1),sen=[sen '*'];end
                    %if(p<.0001),sen=[sen '*'];end
                    %if(p<.00001),sen=[sen '**'];end
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
                    h=text((ii+jj)*.5,or+d2*.2,sen); set(h,'fontsize',FS+szc,'horizontalalignment','center','verticalalignment','bottom');
                    % h=text((ii+jj)*.5,or+d2*.2,sen);set(h,'fontsize',FS+szc,'horizontalalignment','center','verticalalignment','bottom');
                end
                ocu(ii:jj,1+round((or-oc0)/(ocF*d2)))=1;
                
            end
        end
    end
end



end