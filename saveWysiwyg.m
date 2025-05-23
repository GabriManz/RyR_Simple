function quote = saveWysiwyg( h, filename, paperType )

%SAVEWYSIWYG Save Figure and preserve proportions (What You See Is What You
%Get)
%   saveWysiwyg(h, filename, paperType)
%   Will save the Figure with handle h to file called filename.
%   The format of the file is determined from the extension of FILENAME.
%   
% Author:  William Buller, Michigan Tech Research Institute
% william.buller@mtu.edu
%

%% Default Page Size
% Defeault option, 'Tight', makes the whole page the size of the figure.
if(nargin < 3)
	paperType = 'Tight';
end

p = get(h,'Position');
p(1:2)=0;

set(h,'PaperUnits', 'points')

% Additional options set the figure to the size of US_Letter, or A4 in
% Portrait Orientation.  
switch paperType
	case 'Tight'
		pageWidth = p(3);
		pageHeight = p(4);
		marginLeft = 0;
		marginTop = 0;
	case 'US_Letter'
		pageWidth = 612;
		pageHeight = 792;
		marginLeft = 72;
		marginTop = 72;
		p = fitplot(p,pageWidth,pageHeight,marginLeft,marginTop);
	case 'A4'
		pageWidth = 595;
		pageHeight = 842;
		marginLeft = 57;
		marginTop = 57;
		p = fitplot(p,pageWidth,pageHeight,marginLeft,marginTop);
end

set(h, 'PaperPosition', p);
saveas(h,filename);

quote = 'Commonsense is the realised sense of proportion. -Mohandas Gandhi';

end

function p = fitplot(p,pageWidth,pageHeight,marginLeft,marginTop)		
maxWidth = pageWidth-2*marginLeft;
maxHeight = pageHeight-2*marginTop;
p(3:4) = p(3:4)*maxWidth/p(3);
if(p(4)>maxHeight)
p(3:4) = p(3:4)*maxHeight/p(4);
end

end