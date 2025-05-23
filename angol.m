function [a]=angol(v1,v2)
% retorna l'angle positiu (antihorari) entre el vector primer i el segon

x0=v1(1);y0=v1(2);
xf=v2(1);yf=v2(2);

a=acos((x0*xf+y0*yf)/norm([x0 y0])/norm([xf yf]));

a=a*180/pi;

c=cross([x0 y0 0],[xf yf 0]);
c=c(3);

if(c<0), a=360-a;end

end


