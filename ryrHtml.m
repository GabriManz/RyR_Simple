function [] = ryrHtml(folder, RF, Q,Qr,th, fr, tn, roiR, thI, minR, maxR, as)

ind=find(folder=='\',1,'last');
titol=folder(ind+1:end);


meta={};
meta1={};
cb=0;
cb=cb+1;meta{cb,1}=['Number of images'];
cb=cb+1;meta{cb,1}=['Number of frames'];
cb=cb+1;meta{cb,1}=['Global threshold'];
cb=cb+1;meta{cb,1}=['ROI radius'];
cb=cb+1;meta{cb,1}=['I threshold'];
cb=cb+1;meta{cb,1}=['Min radius'];
cb=cb+1;meta{cb,1}=['Max radius'];
cb=cb+1;meta{cb,1}=['Zline angle'];

cb=0;
cb=cb+1;meta1{cb,1}=[num2str(tn)];
cb=cb+1;meta1{cb,1}=[num2str(fr)];
cb=cb+1;meta1{cb,1}=[num2str(th)];
cb=cb+1;meta1{cb,1}=[num2str(roiR)];
cb=cb+1;meta1{cb,1}=[num2str(thI)];
cb=cb+1;meta1{cb,1}=[num2str(minR)];
cb=cb+1;meta1{cb,1}=[num2str(maxR)];
cb=cb+1;meta1{cb,1}=[num2str(mean(as))];


cc=0;
cc=cc+1;parrafada{cc}='<html>';
cc=cc+1;parrafada{cc}='<head runat="server">';
cc=cc+1;parrafada{cc}='<style>';
cc=cc+1;parrafada{cc}='table {text-align: center; font-size: 12}';
cc=cc+1;parrafada{cc}='h2 {text-align: center; margin: 3px 3px 3px 3px;}';
cc=cc+1;parrafada{cc}='h3 {text-align: center; margin: 3px 3px 3px 3px;}';
cc=cc+1;parrafada{cc}='h4 {text-align: center; margin: 3px 3px 3px 3px;}';
cc=cc+1;parrafada{cc}='#header {background-color:black;color:white;text-align:center;padding:1px;max-height:4%;}';
cc=cc+1;parrafada{cc}='#head {background-color:black;padding:5px;max-height:6%;float:center;}';
cc=cc+1;parrafada{cc}='#cont0 {overflow-y:auto;border: 3px solid #000000;background-color:#eeeeee;height:90%;width:40%;float:left;}';
cc=cc+1;parrafada{cc}='#cont1 {overflow-y:auto;border: 3px solid #000000;background-color:#eeeeee;height:90%;width:58%;float:right;}';
cc=cc+1;parrafada{cc}='#nav1 {overflow-y:auto;padding:5px;background-color:#99ee99;width:98%;float:left;}';
%cc=cc+1;parrafada{cc}='#nav2 {overflow-y:auto;padding:5px;background-color:#ffffaa;max-width:100%;float:center;}';
cc=cc+1;parrafada{cc}='#nav3 {overflow-y:auto;padding:5px;background-color:#ee9999;width:98%;float:right;}';
cc=cc+1;parrafada{cc}='</style>';
cc=cc+1;parrafada{cc}='</head>';

I=pij();imwrite(I,[folder '/' RF '/LASSIE1.png']);
extrashit=['<a href="http://lassielab.cat" ><img src="LASSIE1.png" alt="lassie" style="float:right;height:45px;padding: 10px;"></a>'];
cc=cc+1;parrafada{cc}=['<div id="header"><h2>' titol extrashit '</h2></div>'];%style="width:304px;height:228px;"
cc=cc+1;parrafada{cc}=['<div id="head">'];
%%% taula general data
cc=cc+1;parrafada{cc}=['<table border="0">'];
for ii=1:length(meta)%-1
    cc=cc+1;parrafada{cc}=['<col width="' num2str(floor(100/2/length(meta))) '%">'];
    cc=cc+1;parrafada{cc}=['<col width="' num2str(floor(100/2/length(meta))) '%">'];
    cc=cc+1;parrafada{cc}=['<col width="1">'];
end
cc=cc+1;parrafada{cc}=['<col width="' num2str(floor(100/2/length(meta))) '%">'];
cc=cc+1;parrafada{cc}=['<col width="' num2str(floor(100/2/length(meta))) '%">'];
cc=cc+1;parrafada{cc}='<tr>';
for ii=1:length(meta)
    if(ii==1),
        cc=cc+1;parrafada{cc}=['<th bgcolor=#ffffff>' meta{ii} '</th>']; 
        cc=cc+1;parrafada{cc}=['<th bgcolor=#ffffff>' meta1{ii} '</th>'];
    else
        cc=cc+1;parrafada{cc}=['<th bgcolor=#000000> </th>'];
        cc=cc+1;parrafada{cc}=['<th bgcolor=#ffffff>', meta{ii}, '</th>'];
        switch ii
            case {3,4,5,6,7,8}
                fr=niceNums(str2num(meta1{ii}),2);
            otherwise
                fr=meta1{ii};
        end
        cc=cc+1;parrafada{cc}=['<th bgcolor=#ffffff>' fr '</th>'];
    end
end
cc=cc+1;parrafada{cc}='</tr>';

cc=cc+1;parrafada{cc}='</table>';
cc=cc+1;parrafada{cc}='</div>';


%%%% taules
cc=cc+1;parrafada{cc}=['<div id="cont0">'];

for ii=1:length(Q)
    
    dades=Q{ii};
    dadesOut=Qr{ii};
cc=cc+1;parrafada{cc}=['<h4><center> Image ',num2str(ii),'</center></h5>'];
%%% taula ryrs
cc=cc+1;parrafada{cc}=['<table>'];%<caption>Accepted sparks</caption>'];
cc=cc+1;parrafada{cc}='<tr>';
cc=cc+1;parrafada{cc}='<th valign="top">';
cc=cc+1;parrafada{cc}=['<div id="nav1">'];
cc=cc+1;parrafada{cc}=['<table border="1">'];

cc=cc+1;parrafada{cc}='<tr>';
for ii=1:size(dades,2)
    cc=cc+1;parrafada{cc}=['<col width="80">'];
end
cc=cc+1;parrafada{cc}='<tr>';
cc=cc+1;parrafada{cc}=['<th>RyR ID</th><th>Xpix</th><th>Ypix</th><th>Mean I</th><th>Max I</th><th>R (um)</th><th>DistNN (um)</th><th>Zline</th><th>Layer</th>'];
cc=cc+1;parrafada{cc}='</tr>';
for ii=1:size(dades,1)
    cc=cc+1;parrafada{cc}='<tr>';
    frase=[];
    for jj=1:size(dades,2)
             switch jj
                case {4,5,6,7}
                    fr=niceNums(dades(ii,jj),2);
                otherwise
                    fr=num2str(dades(ii,jj));
            end
         frase=[frase '<td>' fr '</td>'];
    end
    cc=cc+1;parrafada{cc}=frase;
    cc=cc+1;parrafada{cc}='</tr>';
end
%cc=cc+1;parrafada{cc}=['<tr><th>.</th></tr>'];


cc=cc+1;parrafada{cc}='</table>';
cc=cc+1;parrafada{cc}='</div>';

cc=cc+1;parrafada{cc}='</th>';




%%% taula Rejected

cc=cc+1;parrafada{cc}='<th valign="top">';
cc=cc+1;parrafada{cc}=['<div id="nav3">'];
cc=cc+1;parrafada{cc}=['<table border="1">'];


for ii=1:size(dades,2)
    cc=cc+1;parrafada{cc}=['<col width="80">'];
end
cc=cc+1;parrafada{cc}='<tr>';
cc=cc+1;parrafada{cc}=['<th>RyR ID</th><th>Xpix</th><th>Ypix</th><th>Mean I</th><th>Max I</th><th>R (um)</th><th>DistNN (um)</th><th>Zline</th><th>Layer</th>'];
cc=cc+1;parrafada{cc}='</tr>';
for ii=1:size(dadesOut,1)
    cc=cc+1;parrafada{cc}='<tr>';
    frase=[];
    for jj=1:size(dadesOut,2)
            switch jj
                case {4,5,6,7}
                    fr=niceNums(dadesOut(ii,jj),2);
                otherwise
                    fr=num2str(dadesOut(ii,jj));
                    
            end
        frase=[frase '<td>' fr '</td>'];
    end
    cc=cc+1;parrafada{cc}=frase;
    cc=cc+1;parrafada{cc}='</tr>';
end


cc=cc+1;parrafada{cc}=['</table>'];
cc=cc+1;parrafada{cc}='</div>';
cc=cc+1;parrafada{cc}='</th>';

cc=cc+1;parrafada{cc}='</tr>';

cc=cc+1;parrafada{cc}='</table>';

end

cc=cc+1;parrafada{cc}='</div>';



%%% contingut

cc=cc+1;parrafada{cc}=['<div id="cont1"><center>' ];
cc=cc+1;parrafada{cc}='<h3>Cluster Locations</h3>';
for ii=1:length(Q)


    filename=['RyRsImage',num2str(ii),'_',titol];
    cc=cc+1;parrafada{cc}=['<h4> Image ',num2str(ii),'</h4>'];
cc=cc+1;parrafada{cc}=['<div><img id="fotards" src="',filename,'.png" /></div>'];
 cc=cc+1;parrafada{cc}=['<h4> </h4>'];
filename=['RyRsPealing',num2str(ii),'_',titol];
cc=cc+1;parrafada{cc}=['<div><img id="fotards" src="',filename,'.png" /></div>'];
end

cc=cc+1;parrafada{cc}='</center></div>';





cc=cc+1;parrafada{cc}='</html>';

%%%munta html
rutadades=[folder,'\', RF '\' titol '.html'];
pumpum=fopen(rutadades,'wt');
for ii=1:length(parrafada)
    fprintf(pumpum,'%s\n',parrafada{ii});
end
fclose(pumpum);



end
