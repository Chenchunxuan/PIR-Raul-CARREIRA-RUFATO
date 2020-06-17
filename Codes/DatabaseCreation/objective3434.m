function E = objective3434(A,~)

b8=A(1);b0=A(2);b2=A(3);b15=A(4);b17=A(5);Bte=A(6);Yle=A(7);Ate=A(8); 
% a1=A(1);a2=A(2);a3=A(3);a4=A(4);a5=A(5);a6=A(6);a7=A(7);a8=A(8);a9=A(9);
% a10=A(10);a11=A(11);a12=A(12);a13=A(13);a14=A(14);a15=A(15);a16=A(16);
% global X Y yt xt yc xc dzte zte rle m n
 [X,Y,yt,xt,yc,xc,dzte,zte,rle,m,n] = fcn_2();
%% Bezier polynomials

[LET,LEC,TET,TEC] = CurvePoints(xt,yt,xc,yc,Bte,Yle,Ate,dzte,zte,rle,b8,b0,b2,b15,b17);
% LET=[0 0 a1 xt;0 a2 yt yt];
% LEC=[0 a3 a4 xc;0 a5 yc yc];
% TET=[xt a6 a7 a8 1;yt yt a9 a10 0];
% TEC=[xc a12 a13 a14 1;yc yc a15 a16 0];

Plet=[-1*LET(:,1)+3*LET(:,2)-3*LET(:,3)+1*LET(:,4) 3*LET(:,1)+-6 ...
    *LET(:,2)+3*LET(:,3) -3*LET(:,1)+3*LET(:,2) LET(:,1)];

Plec=[-1*LEC(:,1)+3*LEC(:,2)-3*LEC(:,3)+1*LEC(:,4) 3*LEC(:,1)+-6 ...
    *LEC(:,2)+3*LEC(:,3) -3*LEC(:,1)+3*LEC(:,2) LEC(:,1)];

Ptet=[1*TET(:,1)-4*TET(:,2)+6*TET(:,3)-4*TET(:,4)+1*TET(:,5) -4 ...
    *TET(:,1)+12*TET(:,2)-12*TET(:,3)+4*TET(:,4) 6*TET(:,1)-12 ...
    *TET(:,2)+6*TET(:,3) -4*TET(:,1)+4*TET(:,2) 1*TET(:,1)];

Ptec=[1*TEC(:,1)-4*TEC(:,2)+6*TEC(:,3)-4*TEC(:,4)+1*TEC(:,5) -4 ...
    *TEC(:,1)+12*TEC(:,2)-12*TEC(:,3)+4*TEC(:,4) 6*TEC(:,1)-12 ...
    *TEC(:,2)+6*TEC(:,3) -4*TEC(:,1)+4*TEC(:,2) 1*TEC(:,1)];

PXlet = Plet(1,:);PYlet=Plet(2,:);
PXlec = Plec(1,:);PYlec=Plec(2,:);
PXtet = Ptet(1,:);PYtet=Ptet(2,:);
PXtec = Ptec(1,:);PYtec=Ptec(2,:);

%% Polynomial points + error function

E  = 0;   % Error
T  = zeros(n,1);
C  = zeros(n,1);
Yn = zeros(n,1);
a  = 0;
i = 0;
for i=1:n
    % thickness
    if X(i)<=xt % leading edge
        A=roots(PXlet-[0 0 0 X(i)]);
        for j=1:3
            if imag(A(j))==0 && (real(A(j))<=1 && real(A(j))>=0)
                a=A(j);
            end
        end
        T(i)= polyval(PYlet,a);
    else % trailing edge

        A = roots(PXtet-[0 0 0 0 X(i)]);

        for j=1:4
            if imag(A(j))==0 && (real(A(j))<=1 && real(A(j))>=0)
                a = A(j);
            end
        end
        T(i) = polyval(PYtet,a);
    end
    
    % camber
    if X(i)<=xc % leading edge
        A = roots(PXlec-[0 0 0 X(i)]);
        for j=1:3
            if imag(A(j))==0 && (real(A(j))<=1 && real(A(j))>=0)
                a=A(j);
            end
        end
        C(i) = polyval(PYlec,a);
    else % trailing edge
        A = roots(PXtec-[0 0 0 0 X(i)]);
        for j=1:4
            if imag(A(j))==0 && (real(A(j))<=1 && real(A(j))>=0)
                a=A(j);
            end
        end         
        C(i) = polyval(PYtec,a);
    end
    if i<=m
        Yn(i)=C(i)+T(i)/2;
    else
        Yn(i)=C(i)-T(i)/2;
    end
    E = E+abs(Yn(i)-Y(i))^2;
end
% figure;plot(X,Y,'o',X,Yn);axis equal

end
