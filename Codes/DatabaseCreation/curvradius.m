function [rle] = curvradius(xa,ya,xb,yb,xc,yc)

m3 = -(xb-xa)/(yb-ya);
m4 = -(xb-xc)/(yb-yc);

n3 = (yb+ya-m3*(xb+xa))/2;
n4 = (yb+yc-m4*(xb+xc))/2;

xc = -(n4-n3)/(m4-m3);
yc = m4*xc+n4;

rle = sqrt((yc-ya)^2+(xc-xa)^2);