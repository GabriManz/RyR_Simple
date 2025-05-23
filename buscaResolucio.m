function [DX,DT,DZ,DY] = buscaResolucio(folder,tagRyR)
% FUNCIONA PER IMATGES I ESCANEJAT DiIMATGES!
% Per DT de LS s'ha d'utilitzar el DY*10^-3 (ms)

namesR=[dir([folder,'\*' tagRyR '*.tif']);dir([folder,'\*' tagRyR '*.png']);dir([folder,'\*' tagRyR '*.jpg']);dir([folder,'\*' tagRyR '*.data'])]; % Imatges RyRs

file=namesR(1).name;

DX=0;DZ=0;DT=0;
try
    % mira de llegir espai i temps del nom darixu
    [DX,DT] = getRes([folder file]);
    
    if(DX==0),
        
        xml=dir([folder '/*.xml']);
        xml=xml(1);
        try, % mira de llegir espai x i temps dun xml
            [DX DY DT TTIME NXPIX NYPIX] = physical_parameters(folder,xml);
        end
        
        % mira de llegir espai x i z dun xml
        [res] = getResolution([folder '/' xml.name]);
        DX=res(1); DZ=res(2);
        if(DX==0),
            
            
            % mira de llegir espai x i temps dun xml
            [DX DY DT TTIME NXPIX NYPIX] = physical_parameters(folder,xml);
            if DT>90, DT=0; end
            if(isempty(DT))
                error('Can''t read physical paremeters');
            end
        end
    end
catch
    error('Can''t read physical paremeters');
end




end