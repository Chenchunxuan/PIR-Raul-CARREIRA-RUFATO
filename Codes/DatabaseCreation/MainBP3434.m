% MAIN PROGRAM FOR AIRFOIL PARAMETRIZATION USING BP3434 CURVES
% Database construction

close all;clear;clc

% Chosing global variables

%% Airfoils characteristics
Nb     = 1000;                     % number of diferent airfoils for tranning data
chord  = 1.0;                     % Airfoils chord
np      = 59;                      % Number of points in airfoil
yt_vec = linspace(0.1,0.4,Nb);   % maximum thickness
yc_vec = linspace(0.05,0.095,Nb); % maximum camber
xc_vec = linspace(0.4,0.9,Nb);  % location of maximmum camber
airfoils = zeros(np,Nb*2);
parameters = zeros(Nb,16);
xt_vec = zeros(1,Nb);

for g=1:Nb
    yt = yt_vec(g);
    yc = yc_vec(g);
    xc = xc_vec(g);
    [coord,xt_vec(g)] = BuildAirfoils(yt,yc,xc,np*0.5+1);
    X   = coord(:,1)';
    Y   = coord(:,2)';
    airfoils(:,2*g-1) = X;
    Xvec(:,g) = X;
    airfoils(:,2*g) = Y;
    Yvec(:,g) = Y;
end
    
NP       = 100; % number of members of the population
D        = 8;   % number of variables to optimize
F        = 0.9; % step between [0 2]
CR       = 0.8; % crossover probability between [0 1]
itermax  = 50;  % maximum number of iterations (generations)
strategy = 7;   % optimization strategy: see options in devec.m
refresh  = 10;  % number of generations
VTR      = 1e-6;% desirable value for objective function


delete(gcp('nocreate'))
parpool(2)
p = 0;

parfor g = 1:Nb
     fprintf('Iteraçao'); disp(g)
    
    % figure;plot(X,Y,'o');axis equal
    % Set the global "VALUE"
    fcn_1(Xvec(:,g),Yvec(:,g),yt_vec(g),xt_vec(g),yc_vec(g),xc_vec(g),np);
    % Read the global "VALUE"
    [X,Y,yt,xt,yc,xc,dzte,zte,rle,m,n] = fcn_2();
    
    %% Definitions for optimization
    
    % b8  = 0 < b8 < min([yt sqrt(-2*rle*xt/3)])
    % b0  = 0 < b0 < xc
    % b2  = b0 < b2 < xc
    % b15 = 3*xt+15*b8^2/(4*rle) < b15 < 1
    % b17 = (-8*yc*cot(Yle)+13*xc)/6 < b17 < 1
    
    LIM = [0 yt;0 xc;0 xc;xt 1;xc 1;0+eps pi/2;0+eps pi/2;0+eps pi/2];
    [bestmem,bestval,nfeval] = devec('objective3434',VTR,D,LIM(:,1)',LIM(:,2)',[],NP,itermax,F,CR,strategy,refresh);
    
    %% Show optimized values
    
    % LET=[0 0 a1 xt;0 a2 yt yt];
    % LEC=[0 a3 a4 xc;0 a5 yc yc];
    % TET=[xt a6 a7 a8 1;yt yt a9 a10 0];
    % TEC=[xc a12 a13 a14 1;yc yc a15 a16 0];
    b8=bestmem(1);b0=bestmem(2);b2=bestmem(3);b15=bestmem(4);b17=bestmem(5);
    Bte=bestmem(6);Yle=bestmem(7);Ate=bestmem(8);
    [LET,LEC,TET,TEC]=CurvePoints(xt,yt,xc,yc,Bte,Yle,Ate,dzte,zte,rle,b8,b0...
        ,b2,b15,b17);
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
    PXlet=Plet(1,:);PYlet=Plet(2,:);
    PXlec=Plec(1,:);PYlec=Plec(2,:);
    PXtet=Ptet(1,:);PYtet=Ptet(2,:);
    PXtec=Ptec(1,:);PYtec=Ptec(2,:);
    E=0;
    T=zeros(n);
    C=zeros(n);
    Yn=zeros(n);
    a=0;
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
    %     1000*E;
%     u=linspace(0,1,100);
%     xC1=zeros(length(u));xT1=xC1;yC1=xC1;yT1=xC1;
%     xC2=xC1;xT2=xC1;yC2=xC1;yT2=xC1;
%     
%     for i=1:100
%         xC1(i)=polyval(Plec(1,:),u(i));xC2(i)=polyval(Ptec(1,:),u(i));
%         xT1(i)=polyval(Plet(1,:),u(i));xT2(i)=polyval(Ptet(1,:),u(i));
%         yC1(i)=polyval(Plec(2,:),u(i));yC2(i)=polyval(Ptec(2,:),u(i));
%         yT1(i)=polyval(Plet(2,:),u(i));yT2(i)=polyval(Ptet(2,:),u(i));
%     end
    
%     figure;plot([xC1 xC2],[yC1 yC2],'-b',[xT1 xT2],[yT1 yT2],'-r');
%     grid on;axis equal;hold on;plot([LEC(1,:),TEC(1,:)],[LEC(2,:),TEC(2,:)],'ob',...
%         [LET(1,:),TET(1,:)],[LET(2,:),TET(2,:)],'or')

%     plot(X,Y,'o',X,Yn);axis equal;legend('Real','Ajusted')
%     drawnow 
    parameters(g,:) = [g xt yt xc yc Bte Yle Ate dzte zte rle b8 b0 b2 b15 b17];
%     plot(X,Y,'o',X,Yn);axis equal;legend('Real','Ajusted')
%     drawnow
    
end

parametersaux = zeros(Nb,15);
for i = 1:Nb

    for j = 1:Nb
        if parameters(j,1)==i
            parametersaux(i,:) = parameters(j,2:end);
        end
    end 

end 

save('Airfoils.txt','airfoils','-ascii')
[fid, msg] = fopen('Parameters.txt', 'w');
for i = 1:length(parametersaux(:,1))
    for  j =1:length(parametersaux(1,:))-1
        s = strcat(num2str(parametersaux(i,j),'%12.16f'),';');
        fprintf(fid, strrep(s, '%12.16f', 'char'));
    end
    fprintf(fid, strrep(num2str(parametersaux(i,length(parametersaux(1,:))),'%12.16f'), '%12.16f', 'char'));
    fprintf(fid,'%c\n',' ');
end
fclose(fid);



