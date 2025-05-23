function [rutes,N] = inicialitza(N,folder)

if(nargin>1),
    
    rutes{1}=folder;
    N=1;
    
else
    if(N>0)
        folder=' ';
        for ii=1:N
            
            folder=uigetdir(folder);
            
            rutes{ii}=folder;
        end
    else
      
            folder=' ';
            folder=uigetdir(folder);
       
        s=dir(folder);
        cc=0;
        for ii=3:length(s)
            if(s(ii).isdir)
                cc=cc+1;
                rutes{cc}=[folder '/' s(ii).name];
            end
        end
        N=length(rutes);
        
    end
end
save('ultimPaquet','rutes');


end