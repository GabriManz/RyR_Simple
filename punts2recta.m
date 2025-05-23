function [m,n] = punts2recta(x1,y1,x2,y2)


m=(y1-y2)/(x1-x2);
n=y1-m*x1;

end