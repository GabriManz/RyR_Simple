function [Z] = escriuFrase(frase,tam)
if(nargin<2),tam=16;end

    
   Z=[];
for ii=1:length(frase)
    
    ruut=['charSet/' frase(ii) '.png'];
    if(exist(ruut)~=2)
       aux='0';
       if(frase(ii)==lower(frase(ii))),aux='2';end
       if(frase(ii)==upper(frase(ii))),aux='1';end
       ruut=['charSet/' frase(ii) aux '.png'];
    
    end
    if(frase(ii)=='µ')
        ruut=['charSet/10.png'];
    end    
    if(exist(ruut)==2)
  X=imread(ruut);
  Z=[Z,zeros(size(X,1),1),X];
    else
        X=zeros(22,6);
   if(frase(ii)=='.')
       X(16:17,2:3)=1;
   end 
   if(frase(ii)==',')
       X(16:17,2:3)=1;
       X(18,2)=1;X(19,2)=1;X(19,1)=1;
   end 
   if(frase(ii)=='-')
       X(11:12,1:end-1)=1;
   end    
   if(frase(ii)=='_')
       X(16:17,1:end-1)=1;
   end    
   if(frase(ii)=='+')
       X(11:12,1:6)=1;
       X(9:14,3:4)=1;
   end
   if(frase(ii)=='=')
       X(9:10,1:6)=1;
       X(13:14,1:6)=1;       
   end
   if(frase(ii)=='(')
       X=imread(['charSet/pL.png']);
   end 
   if(frase(ii)==')')
       X=imread(['charSet/pR.png']);
   end 
   if(frase(ii)=='ñ')
       X=imread(['charSet/n2.png']);
       X(4,2)=1;X(4,5:6)=1;X(3,7)=1;X(3,3:4)=1;
   end 
   if(frase(ii)=='%')
       X=imread(['charSet/x2.png']);
      X(8:9,2:4)=0;X(10:11,3:4)=0;
      X(13:15,6:8)=0;X(7,3)=0;
      X(6:11,:)=X(4:9,:);
      X(1:3,:)=0;
   end    
   if(frase(ii)=='!')
       X=imread(['charSet/i2.png']);
       X(4:5,1)=1;X(14:15,1)=0;
   end 
   if(frase(ii)=='?')
       X=imread(['charSet/qM.png']);
   end   
   if(frase(ii)=='<')
       X=imread(['charSet/lT.png']);
   end   
   if(frase(ii)=='>')
       X=imread(['charSet/gT.png']);
   end  
   if(frase(ii)=='¿') % interrogant invertit vol dir 'a la menys 1'
X(4,1:2)=1;X(2,3)=1;X(2:7,4)=1;
   end   
   if(frase(ii)=='²')
X(3,1)=1;X(2,2)=1;X(3:4,3)=1;X(5,2)=1;X(6:7,1)=1;X(7,2:3)=1;     
   end
   if(frase(ii)=='³')
X(2,1)=1;X(2,2)=1;X(3,3)=1;X(4,2)=1;X(5:6,3)=1;X(7,2)=1;X(6,1)=1;
   end
   Z=[Z,zeros(size(X,1),1),X];
    end
end



if(tam~=16),
    if(tam<12)
    Zs=bwmorph(Z,'skel',Inf);
    [a,b]=find(Zs==1);
    a=floor(a*tam/16);
    b=floor(b*tam/16);
    out=[find(a==0) ;find(b==0)];
    a(out)=[];b(out)=[];
    Z1=zeros(size(Z));
    for ii=1:length(a)
        Z1(a(ii),b(ii))=1;
    end
    Z=Z1(1:ceil(size(Z1,1)*tam/16),1:ceil(size(Z1,2)*tam/16));
   % figure;imagesc(Z1),axis equal
    else
    Z=imresize(Z,tam/16);
    end
end

% [a,b]=find(Z==1);
% 
%     
% mr=round(tam/7);if(mr<1),mr=1;end
% if(sum(a<=mr)>0),Z=[zeros(mr,size(Z,2));Z];end
% if(sum(b<=mr)>0),Z=[zeros(size(Z,1),mr),Z];end
% [a,b]=find(Z==1);
% if(sum(a+2*mr>size(Z,1))>0),Z=[Z;zeros(2*mr,size(Z,2))];end
% if(sum(b+2*mr>size(Z,2))>0),Z=[Z,zeros(size(Z,1),2*mr)];end
% 
% Z=Z(min(a)-mr:max(a)+2*mr,min(b)-mr:max(b)+2*mr);

end






% 
% 
% frase='0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
% 
% baix='gpqy';
% %dalt='bdhiklt'
% %tots2='fj';
% frase=baix;
% for ii=1:length(frase)
%     
%  ruut=['charSet/' frase(ii) '.png'];
%     if(exist(ruut)~=2)
%        aux='0';
%        if(frase(ii)==lower(frase(ii))),aux='2';end
%        if(frase(ii)==upper(frase(ii))),aux='1';end
%        ruut=['charSet/' frase(ii) aux '.png'];
%     
%     end
%  X=imread(ruut);
%  X(end-1:end,:)=[];X=[zeros(2,size(X,2));X];imwrite(X,ruut);
%  %if(size(X,1)==20),X=[X;zeros(2,size(X,2))];imwrite(X,ruut);end
% disp(size(X));
% end
% 


