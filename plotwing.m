% H = plotwing(BLADE,XSECTION) ...
% H = plotwing(BLADE,XSECTION,AOA) ...
% 
% Generates a 3D plot of a wing/blade. BLADE contains the filename where
% the geometric information is stored.  XSECTION contains the shape of the
% airfoil profile.  This function cannot handle blade designs that
% transition between two airfoils.
% 
% NOTE: This function expects the twist angle in the BLADE file to be
% expressed in degrees.  The data should be stored in columns, giving
% radial position (r), chord (c) and twist (theta), in that order.
%
% AOA is an optional input that sets a collective angle of attack for the
% whole wing, expressed in degrees, nose-up positive.
% 
%  ---------------------------------------------------------------------

function h = plotwing(blade,xsection,varargin)

    if ~isempty(varargin)
        alfa = varargin{1};
    else 
        alfa = 0;
    end

    RiCi = readmatrix(blade);
    XiYi = readmatrix(xsection);

    Ri = RiCi(:,1);
    Ci = RiCi(:,2);
    Bi = RiCi(:,3)*pi/180; %********************************************
    %----------------------------------------------------------------
        nx = length(XiYi(:,1));         nc = length(Ci);
        ny = length(XiYi(:,2));         nr = length(Ri);
    %----------------------------------------------------------------
        minr = min(Ri);
        maxr = max(Ri);
    %...............................................................
        [xm,ym] = perfil(xsection);

    %...............................................................
    %...............................................................
    % INTERPOLA AS COORDENADAS DOS PONTOS DO PERFIL EM CADA ESTA�AO
    %---------------------------------------------------------------
       xmin = min(XiYi(:,1)-xm);       ymin = min(XiYi(:,2)-ym);
       xmax = max(XiYi(:,1)-xm);       ymax = max(XiYi(:,2)-ym);

       dx0 = (xmax-xmin)/(nx-1);       x0 = xmin:dx0:xmax;
       dy0 = (ymax-ymin)/(ny-1);       y0 = ymin:dy0:ymax;

       dx = (xmax-xmin)/(nr-1);        x = xmin:dx:xmax;
       dy = (ymax-ymin)/(nr-1);        y = ymin:dy:ymax;

       xx = spline(x0,XiYi(:,1)-xm,x);
       yy = spline(y0,XiYi(:,2)-ym,y);

    %...............................................................
    % ESTABELECE AS COORDENADAS DOS PONTOS DO PERFIL EM CADA ESTA�AO
    %---------------------------------------------------------------
    %             FX = maxr*Ci*xx;
    %			 FY = maxr*Ci*yy;

                 FX = Ci*xx;
                 FY = Ci*yy;

    %.............................................................
    % PROVOCA A ROTA�AO DOS VETORES QUE FORMAM AS ESTA�OES DA PA
    %-------------------------------------------------------------
            for i = 1:nr
                FXr(:,i) = cos(Bi).*FX(:,i) - sin(Bi).*FY(:,i);
                FYr(:,i) = sin(Bi).*FX(:,i) + cos(Bi).*FY(:,i);
            end

    %--------------------------------------------------------------------------
    %for i=1:nr
    %    plot(FX(i,:),FY(i,:));grid on;axis equal;hold on
    %end
    %return

    z = linspace(minr,maxr,nr);

    r = sqrt(FXr.*FXr + FYr.*FYr);

    theta = atan2(FYr,FXr);

        Ax = r.*cos(theta-alfa*pi/180);
        Ay = r.*sin(theta-alfa*pi/180);
        Az = z'*ones(1,nr);

    h = figure();
    mesh(Ax,Az,Ay); axis equal;
    set(gca,'FontName','times','FontSize',14)
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis on; grid off;  box on;
    colormap([0 0 0]);
    view(-45,45);
%     view(90,45); axis off
    rotate3d on;

end