function [] = ensenya(frase)


   

cl=clock;
if(cl(4)<10) cl4=['0' num2str(cl(4))];else cl4=num2str(cl(4)); end 
if(cl(5)<10) cl5=['0' num2str(cl(5))];else cl5=num2str(cl(5)); end 
cl6=floor(cl(6)); 
if(cl6<10) cl6=['0' num2str(cl6)];else cl6=num2str(cl6); end 
if(frase(1)~=' ')
disp([cl4 ':' cl5 ':' cl6 ' ' frase]);    
else
disp([' ' cl4 ':' cl5 ':' cl6 ' ' frase]); 
end    
end