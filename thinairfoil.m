% [CL,CM] = THINAIRFOIL(AIRFOILNAME,ALFA) ...
% 
% Lift and moment coefficients (about the quarter-chord) are calculate for
% a thin cambered airfoil.  
% 
% AIRFOILNAME is the filename where the airfoil geometry is stored in Selig
% format.  This is NOT the camber line, but rather is the geometry
% including points on the pressure and suction sides.  The camber line is
% calculated from these data.
% 
% ALFA is the horizontal vector of angles of attack at which CL and CM are
% to be calculated, expressed in degrees.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [cl,cm] = thinairfoil(airfoilname,alfa)

    % Load airfoil geometry file
    data = readmatrix(airfoilname);
    x = data(:,1);
    y = data(:,2);
    
    % Split airfoil coordinates into upper and lower surface
    indle = find(x==min(abs(x))); % index of leading-edge coordinates in x vector
    xu = flipud(x(1:indle));
    yu = flipud(y(1:indle));
    xl = x(indle:end);
    yl = y(indle:end);

    % Calculate camber line
    xc = xu;
    ycl = interp1(xl,yl,xu);
    yc = 0.5*(yu + ycl);

    % Calculate gamma distribution for each angle of attack
    alfa = pi/180*alfa; % angle of attack, converted to radians
    As = fouriercoeffs(xc,yc,alfa);
    gmmas = zeros(length(xc),3);
    for ii = 1:length(alfa)
        gmmas(:,ii) = gammacalc(xc,yc,As(:,ii));   
    end

    % Calculate lift coefficient
    cl = zeros(size(alfa));
    for ii = 1:length(alfa)
        cl(ii) = 2*pi*(As(1,ii) + As(2,ii)/2);
    end

    % Calculate moment coefficient
    cm = zeros(size(alfa));
    for ii = 1:length(alfa)
        cm(ii) = pi/4 * (As(3,ii) - As(2,ii));
    end


end


%%%%%%

function A = fouriercoeffs(xc,yc,alfa)

    % Calculate the transformation variable, thta
    thta = theta(xc);

    % Calculate dz/dx 
    dzdx = gradient(yc)./gradient(xc);

    % Initialize A matrix
    a = round(length(xc)); 
    b = length(alfa);
    A = zeros(a,b);

    % leading coefficient, A0 in notes
    A(1,:) = alfa - (1/pi) * ones(size(alfa)) * trapz(thta, dzdx);

    % other Fourier coefficients of vorticity distribution
    for ii = 2:a
        n = ii-1;
        A(ii,:) = (2/pi) * ones(size(alfa))* trapz(thta,dzdx.*cos(n*thta));
    end
    
end


function gmma = gammacalc(xc,yc,A)
% Calculate gamma distribution from Fourier coefficients

    % Calculate the transformation variable, thta
    thta = theta(xc);

    % Calculate gamma distribution, eq. (4.47)
    gmma = 2*A(1).*(1+cos(thta))./sin(thta);
    for ii = 2:length(A)
        n = ii - 1;
        gmma = gmma + 2*A(ii) * sin(n*thta);
    end

end

function thta = theta(xc)

    thta = acos(1-2*xc);

end