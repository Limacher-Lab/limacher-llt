function h = plotcamberline(airfoilname)

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
    
    % Plot airfoil shape to scale
    h = plot(x,y,'b.',xc,yc,'r.');
    axis equal;
    grid on;
    set(gcf,'color','w');
    legend('airfoil coordinates','camber-line coordinates',...
        'location','northoutside');

end