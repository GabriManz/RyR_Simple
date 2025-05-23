function [ori1,des1,ori2,des2]=pintaCapCul(Cul,Cap,tamany,col,gruix)

 
%%%%%%%%%% CARMETA
unitari=(Cap-Cul)/norm(Cap-Cul);


uP1=[unitari(2),-unitari(1)];
uP2=[-unitari(2),unitari(1)];

% pinta cap
ori=[Cap-unitari*tamany];
des=Cap;
% h1=line([ori(1),des(1)],[ori(2),des(2)]);
%  set(h1,'color',col,'linewidth',gruix);
ori1=des;
des1=[Cap-unitari*tamany+uP1*tamany];
h1=line([ori1(1),des1(1)],[ori1(2),des1(2)]);
 set(h1,'color',col,'linewidth',gruix);
ori=des;
des=[Cap-unitari*tamany+uP2*tamany];
% h1=line([ori(1),des(1)],[ori(2),des(2)]);
%  set(h1,'color',col,'linewidth',gruix);
ori2=des;
des2=Cap;
h1=line([ori2(1),des2(1)],[ori2(2),des2(2)]);
 set(h1,'color',col,'linewidth',gruix);
 
%%%%%%%%%% CARMETA

 
%%%%%%%%%% ALEX
% unitari=(Cap-Cul)/norm(Cap-Cul);
% 
% 
% uP1=[unitari(2),-unitari(1)];
% uP2=[-unitari(2),unitari(1)];
% 
% % pinta cap
% ori=[Cap-unitari*tamany];
% des=Cap;
% % h1=line([ori(1),des(1)],[ori(2),des(2)]);
% %  set(h1,'color',col,'linewidth',gruix);
% ori=des;
% des=[Cap-unitari*tamany+uP1*tamany];
% h1=line([ori(1),des(1)],[ori(2),des(2)]);
%  set(h1,'color',col,'linewidth',gruix);
% ori=des;
% des=[Cap-unitari*tamany+uP2*tamany];
% % h1=line([ori(1),des(1)],[ori(2),des(2)]);
% %  set(h1,'color',col,'linewidth',gruix);
% ori=des;
% des=Cap;
% h1=line([ori(1),des(1)],[ori(2),des(2)]);
%  set(h1,'color',col,'linewidth',gruix);
%%%%%%%%%% ALEX
 
% pinta cul
% Cap=Cul;
% ori=[Cap-unitari*tamany];
% des=Cap;
% h1=line([ori(1),des(1)],[ori(2),des(2)]);
%  set(h1,'color',col,'linewidth',gruix);
% ori=des;
% des=[Cap-unitari*tamany+uP1*tamany];
% h1=line([ori(1),des(1)],[ori(2),des(2)]);
%  set(h1,'color',col,'linewidth',gruix);
% ori=des;
% des=[Cap-unitari*tamany+uP2*tamany];
% h1=line([ori(1),des(1)],[ori(2),des(2)]);
%  set(h1,'color',col,'linewidth',gruix);
% ori=des;
% des=Cap;
% h1=line([ori(1),des(1)],[ori(2),des(2)]);
%  set(h1,'color',col,'linewidth',gruix);


end