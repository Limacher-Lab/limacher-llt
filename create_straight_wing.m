% Create elliptic blade geometry file for aspect ratio AR with N elements.

function create_straight_wing(filename,N,AR)

    y = linspace(-0.5,0.5,N);
    x = (1/AR) * ones(size(y));
    
    th = zeros(size(y));
    
    data = [y',x',th'];
    
    T = array2table(data,'VariableNames',{'y','c','th'});
    writetable(T,filename);

end