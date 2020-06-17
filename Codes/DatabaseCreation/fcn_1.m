function fcn_1(a,b,c,d,e,f,g)
global Xi Yi yti xti yci xci dztei ztei rlei mi ni
Xi = a;
Yi = b;
yti = c;
xti = d;
yci = e;
xci = f;
dztei = Yi(1)-Yi(end); % final thickness
[~, mi] = min(Xi); % Division point Upper surface/Lower Surface
ztei = Yi(1)+Yi(end)/2;  % final height deviatio
rlei = -curvradius(Xi(mi-1),Yi(mi-1),Xi(mi),Yi(mi),Xi(mi+1),Yi(mi+1)); % radius curvature in leading edge

ni = g;
end
