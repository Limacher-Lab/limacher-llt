% Create elliptic blade geometry file for aspect ratio AR with N elements.

function create_elliptic_blade(filename,N,AR)

    y = linspace(-0.5,0.5,N);
    
    a = 4/(pi * AR);

    c = 2*a * (1/4 - y.^2).^0.5;
    
    th = zeros(size(y));
    
    data = [y',c',th'];
    
    T = array2table(data,'VariableNames',{'y','c','th'});
    writetable(T,filename);

end