function [resolution1,resolution2] = getRes(FileName)
%retorna r1 en um i r2 en ms
if(FileName(end)=='\'),FileName=FileName(1:end-1);end
% if(~isempty(find(FileName=='\'))),
% FileName=FileName(find(FileName=='\',1,'last')+1:end);
% end
noms={'first','second','third','fourth','fifth','sixth'};
FileName=lower(['sisvuit' FileName]);

resolution=0;

% search um
conta=0;
aux=strfind(FileName,'um');
if(~isempty(aux))
    
    for ii=1:length(aux)
        fin=aux(ii)-1;
        faf=aux(ii);
        ori=0;
        while(faf>1)
            
            faf=faf-1;
            if(isempty(strfind(' 0123456789',FileName(faf))))
                if(FileName(faf)~='.')
                    ori=faf+1;faf=0;
                    break;
                end
            end
            
            
        end
        if(ori<=fin)
            if(FileName(ori)=='.'),ori=ori+1;end
            conta=conta+1;
            resolution(conta)=str2num(FileName(ori:fin));
        end
    end
    
end



% search micron
aux=strfind(FileName,'micron');
if(~isempty(aux))
    
    for ii=1:length(aux)
        fin=aux(ii)-1;
        faf=aux(ii);
        ori=0;
        while(faf>1)
            
            faf=faf-1;
            if(isempty(strfind(' 0123456789',FileName(faf))))
                if(FileName(faf)~='.')
                    ori=faf+1;faf=0;
                    break;
                end
            end
            
            
        end
        if(ori<=fin)
            conta=conta+1;
            resolution(conta)=str2num(FileName(ori:fin));
        end
    end
    
end


% search nicron
aux=strfind(FileName,'nicron');
if(~isempty(aux))
    
    for ii=1:length(aux)
        fin=aux(ii)-1;
        faf=aux(ii);
        ori=0;
        while(faf>1)
            
            faf=faf-1;
            if(isempty(strfind(' 0123456789',FileName(faf))))
                if(FileName(faf)~='.')
                    ori=faf+1;faf=0;
                    break;
                end
            end
            
            
        end
        if(ori<=fin)
            conta=conta+1;
            resolution(conta)=str2num(FileName(ori:fin));
        end
    end
    
end

% search pix
aux=strfind(FileName,'pix');
if(~isempty(aux))
    
    for ii=1:length(aux)
        fin=aux(ii)-1;
        faf=aux(ii);
        ori=0;
        while(faf>1)
            
            faf=faf-1;
            if(isempty(strfind(' 0123456789',FileName(faf))))
                if(FileName(faf)~='.')
                    ori=faf+1;faf=0;
                    break;
                end
            end
            
            
        end
        if(ori<=fin)
            conta=conta+1;
            resolution(conta)=str2num(FileName(ori:fin));
        end
    end
    
end

if(length(resolution)>1),
cprintf('err',    ['Warning: Found multiple sizes.\n']);
resolution=resolution(1);
end
resolution1=resolution;


resolution=0;

% search ms
conta=0;
aux=strfind(FileName,'ms');
if(~isempty(aux))
    
    for ii=1:length(aux)
        fin=aux(ii)-1;
        faf=aux(ii);
        ori=0;
        while(faf>1)
            
            faf=faf-1;
            if(isempty(strfind(' 0123456789',FileName(faf))))
                if(FileName(faf)~='.')
                    ori=faf+1;faf=0;
                    break;
                end
            end
            
            
        end
        if(ori<=fin)
            if(FileName(ori)=='.'),ori=ori+1;end
            conta=conta+1;
            resolution(conta)=str2num(FileName(ori:fin));
        end
    end
    
end
aux=strfind(FileName,'fps');
if(~isempty(aux))
    
    for ii=1:length(aux)
        fin=aux(ii)-1;
        faf=aux(ii);
        ori=0;
        while(faf>1)
            
            faf=faf-1;
            if(isempty(strfind(' 0123456789',FileName(faf))))
                if(FileName(faf)~='.')
                    ori=faf+1;faf=0;
                    break;
                end
            end
            
            
        end
        if(ori<=fin)
            if(FileName(ori)=='.'),ori=ori+1;end
            conta=conta+1;
            resolution(conta)=1000/str2num(FileName(ori:fin));
        end
    end
    
end
if(length(resolution)>1),
cprintf('err',    ['Warning: Found multiple times.\n']);
resolution=resolution(1);
end
resolution2=resolution;
end