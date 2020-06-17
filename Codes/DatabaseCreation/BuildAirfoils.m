function [coord,xt] = BuildAirfoils(t,m,p,NeL)
%NACA4: Outputs NACA4 Coordinates
% yc: y due to camber
% yt: y due to thickness
% c: chord
% t: maximum thickness
% m: maximum camber
% p: location of maximmum camber
x = linspace(0,1,NeL).^2;
ytt = x*0.0; ycc = x*0.0;

for i = 1:length(x)
    ytt(i) = 5*t*((0.2969*sqrt(x(i))-0.126*(x(i))-0.3516*(x(i))^2+0.2843*(x(i))^3)+(-0.1015*(x(i))^4));
   
    if x(i)<p
        ycc(i) = (m/p^2)*(2*p*x(i)-x(i)^2);
    else
        ycc(i) = (m/((1-p)^2))*((1-2*p)+2*p*x(i)-x(i)^2);
    end
   
    %create points surrounding shape
    coord = zeros(2*length(x)-1,2);
    coord(1:length(x),1)= flip(x);
    y2 = ycc + ytt;
    
    coord(1:length(x),2)= flip(y2);
    coord(length(x):end,1) = x(1:length(x));
    y2 = ycc - ytt;
    
    coord(length(x):end,2) = y2(1:length(x));
end
[M,I] = max(ytt);  % location of maximmum thickness
xt = x(I);
end
