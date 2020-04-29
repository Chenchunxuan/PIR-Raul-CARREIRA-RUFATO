function [LET,LEC,TET,TEC]=CurvePoints(xt,yt,xc,yc,Bte,Yle,Ate,dzte,zte,rle...
    ,b8,b0,b2,b15,b17)
%% Calculation of bezier curve points 3434

% Leading edge thickness curve (LET)
% LET : Matrix 4 rows (points) by 2 columns (x, y) of bezier points
LET(1,3) = -3*b8^2/(2*rle);
LET(1,4) = xt;
LET(2,2) = b8;
LET(2,3) = yt;
LET(2,4) = yt;

% Leading edge camber curve (LEC)
% LEC: Matrix 4 rows (points) by 2 columns (x, y) of bezier points
LEC(1,2) = b0;
LEC(1,3) = b2;
LEC(1,4) = xc;
LEC(2,2) = b0*tan(Yle);
LEC(2,3) = yc;
LEC(2,4) = yc;

% Trailing edge thickness curve (TET)
% TET: Matrix 5 rows (points) by 2 columns (x, y) of bezier points
TET(1,1) = xt;
TET(1,2) = (7*xt+9*b8^2/(2*rle))/4;
TET(1,3) = (3*xt-3*LET(1,3))/4;
TET(1,4) = b15;
TET(1,5) = 1;
TET(2,1) = yt;
TET(2,2) = yt;
TET(2,3) = (yt+b8)/2;
TET(2,4) = dzte+(1-b15)*tan(Bte);
TET(2,5) = dzte;

% Trailing edge camber curve (TEC)
% TEC: Matrix 5 rows (points) by 2 columns (x, y) of bezier points
TEC(1,1) = xc;
TEC(1,2) = (3*xc-yc*cot(Yle))/2;
TEC(1,3) = (-8*yc*cot(Yle)+13*xc)/6;
TEC(1,4) = b17;
TEC(1,5) = 1;
TEC(2,1) = yc;
TEC(2,2) = yc;
TEC(2,3) = 5*yc/6;
TEC(2,4) = zte+(1-b17)*tan(Ate);
TEC(2,5) = zte;
end