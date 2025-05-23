function [DX DY DT TTIME NXPIX NYPIX] = physical_parameters(folder,xmlfile)

%fid = fopen(xmlfile.name); 
%C = textscan(fid, '%s',100);
CC = textread([folder,'/',xmlfile(end).name], '%s',200000);

ss = strfind(CC,'NumberOfElements');
dds = [];
for i=1:length(CC)
    if ss{i} == 1 
    dds = [dds i];
    end;
end;
    
ss2 = strfind(CC,'Length');
dds2 = [];
for i=1:length(CC)
    if ss2{i} == 1 
    dds2 = [dds2 i];
    end;
end;

ss3 = strfind(CC,'Unit');
dds3 = [];
for i=1:length(CC)
    if ss3{i} == 1 
    dds3 = [dds3 i];
    end;
end;

% PIXELS IN X DIRECTION:
xind = find(CC{dds(1)}=='"');
NXPIX = str2num(CC{dds(1)}(xind(1)+1:xind(2)-1));
% Length in X direction (in meters)
xind = find(CC{dds2(1)}=='"');
LX = str2num(CC{dds2(1)}(xind(1)+1:xind(2)-1));
% Pixel size in X direction (in microns)
DX = (LX/NXPIX)*1e6;


% PIXELS IN Y DIRECTION:
xind = find(CC{dds(2)}=='"');
NYPIX = str2num(CC{dds(2)}(xind(1)+1:xind(2)-1));
% Length in X direction (in meters)
xind = find(CC{dds2(2)}=='"');
LY = str2num(CC{dds2(2)}(xind(1)+1:xind(2)-1));
% Pixel size in X direction (in microns)
DY = (LY/NYPIX)*1e6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CARMETA
if length(dds3)>3,
    for jj=3:length(dds3)
        ss=CC{dds3(jj)};
        if ~isempty(strfind(ss,'s')), kk=jj; end
    end
else
    kk=3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DT=[];
try,
% Time interval between frames:
xind = find(CC{dds(kk)}=='"');
NFRAMES = str2num(CC{dds(kk)}(xind(1)+1:xind(2)-1));
% Total time (in seconds)
xind = find(CC{dds2(kk)}=='"');
TTIME = str2num(CC{dds2(kk)}(xind(1)+1:xind(2)-1));
% Time between frames (in ms)
DT = (TTIME/NFRAMES)*1e3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CARMETA
if isempty(DT)||DT==0, kk=3;
    % Time interval between frames:
xind = find(CC{dds(kk)}=='"');
NFRAMES = str2num(CC{dds(kk)}(xind(1)+1:xind(2)-1));
% Total time (in seconds)
xind = find(CC{dds2(kk)}=='"');
TTIME = str2num(CC{dds2(kk)}(xind(1)+1:xind(2)-1));
% Time between frames (in ms)
DT = (TTIME/NFRAMES)*1e3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

