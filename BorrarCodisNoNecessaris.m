[fList,pList] = matlab.codetools.requiredFilesAndProducts('RyRactivation_stpau.m');

folder='D:\CODE\Activation Sant Pau';
S=dir('D:\CODE\Activation Sant Pau');

fListtotal=cell(1,size(S,2));
for ii=1:size(S,1)
    fListtotal{ii}=[folder,'\' S(ii).name];
end

ind=setdiff(fList,fListtotal);
ind=ismember(fListtotal,fList);

% ind=1-ind;
for ii=1:length(ind)
    if ind(ii)==0, delete ([S(ii).name])
    end
end