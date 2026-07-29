% [CL,CD,Y,GAMMA,V] = LIFTINGLINE(GEOMFILE,FORCEFILE) ...
% [CL,CD,Y,GAMMA,V] = LIFTINGLINE(GEOMFILE,FORCEFILE,ALFAS) ...
% [CL,CD,Y,GAMMA,V] = LIFTINGLINE(GEOMFILE,FORCEFILE,ALFAS,'relfactor',R)
% 
% Prandtl's Lifting-Line Theory (LLT) solver for a finite wing.
% See Anderson, Fundamentals of Aerodynamics, 5th ed., Ch. 5.
% 
% INPUTS
% ------
% GEOMFILE is the geometry file for the wing, with columns for
% spanwise position and local chord length (both normalized by total span), 
% and local twist angles (in degrees).
% 
% FORCEFILE contains the airfoil lift and drag data, with columns for angle
% of attack (deg), lift coefficient and drag coefficient (in that order).
% Additional columns will be ignored.
% 
% ALFAS is an optional third argument containing the wing angles of attack
% (deg) at which to perform the calculation.  It should be a horizontal 
% vector. If ALFAS is not specified, a wing angle of attack of zero is
% assumed.
% 
% 'relfactor' is an optional name-value pair controlling relaxation (default
% 0.01).  Smaller values improve stability but require more iterations.
%
% OUTPUTS
% -------
% CL and CD are the calculated lift and induced drag coefficients on the 
% full wing.  They will have the same dimensions as ALFAS.
% 
% Y is the spanwise coordinate vector extracted from GEOMFILE,
% normalized by total span.
% 
% GAMMA holds the spanwise circulation distribution for each angle of
% attack in a column pertaining to the position of that angle of attack in
% the ALFAS vector.  Circulation is normalized by freestream velocity and
% total span.
% 
% V holds the spanwise downwash distribution (normalized by freestream 
% velocity) in columns corresponding to each angle of attack in ALFAS, as 
% with GAMMA. 
%
% EXAMPLE
% -------
% [CL, CD] = liftingline('elliptic_AR8.txt', 'thin_airfoil_data.txt', ...
%                        -10:2:10, 'relfactor', 0.01);
% plot(-10:2:10, CL)
%
% The Python equivalent is available at limacher-llt/liftingline.py.
%
% REFERENCES
% ----------
% Anderson, J. D. Jr. (2011). Fundamentals of Aerodynamics, 5th ed.
%   McGraw-Hill.  See Ch. 5 (pp. 449-470) for lifting-line theory,
%   Eq. (5.69) for elliptic-wing lift slope, and Eq. (5.61) for induced drag.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [CL,CD,y,Gamma,v] = liftingline(geomfile,forcefile,varargin)

    % Default values
    AoA = 0;
    relfactor = 0.01;

    % Read optional input of angles of attack
    if nargin > 2
        if ~ischar(varargin{1})
            AoA = varargin{1};
        end
        for ii = 1:(nargin - 2)
            if ischar(varargin{ii})
                switch varargin{ii}
                    case 'relfactor'
                        relfactor = varargin{ii+1};
                end
            end
        end   
    end

    % Ensure AoA is a horizontal vector
    [tempa,tempb] = size(AoA);
    if tempa > tempb
        AoA = transpose(AoA);
    elseif tempa ~= 1 && tempb ~= 1
        disp('AoA cannot be a matrix.');
        return;
    end
    
    % Read force data
        data = dlmread(forcefile, ',', 1, 0);
    alfaref = data(:,1);
    cl = data(:,2);
    cd = data(:,3);

    % Read blade geometry data
        data = dlmread(geomfile, ',', 1, 0);
    y = data(:,1);
    c = data(:,2); 
    th = data(:,3);

    % Ensure chord and span coordinates are normalized by total span
    L = y(end)-y(1); 
    y = y/L;
    c = c/L;

    %%%%%%%%%%%%%%%%%%
    % Check that lift curve is monotonic
    [out, monoInds] = checkmonotonic(cl);
    if out == 0
        disp('Lift coefficient data is not monotonically increasing.')
        disp('  Data will be truncated to its longest monotonic section.')
    end

    % Truncate forcefile to include only monotonic lift section
    alfaref = alfaref(monoInds);
    cl = cl(monoInds);
    cd = cd(monoInds);
    %%%%%%%%%%%%%%%%%%%%%%

    % Intialize solution variables
    temp = zeros(length(y),length(AoA));
    Gamma = temp; 
    v = temp;
    % alfai = temp;
    CL = zeros(size(AoA));
    CD = zeros(size(AoA));

    % Solve
    for ii = 1:length(AoA)
        
        %%%%%%%%%%%%%%%%%
        % [CL(ii),CD(ii),Gamma(:,ii),v(:,ii)] = ...
        %     LLsingle(y,c,th,alfaref,cl,AoA(ii),relfactor);
        [CL(ii),CD(ii),Gamma(:,ii),v(:,ii)] = ...
            LLsingle(y,c,th,alfaref,cl,cd,AoA(ii),relfactor);
        %%%%%%%%%%%%%%%%%

        
    end

end

%%%%%%%%%%
function [CL,CD,Gamma,v] = LLsingle(y,c,th,alfaref,cl,cd,AoA,relfactor)
    
    % Iterate to solve for circulation distribution
    [Gamma,v,alfai,~] = calcgamma(y,c,th,AoA,alfaref,cl,...
        'relfactor',relfactor);

    % Calculate lift and drag
    cli = interp1(alfaref,cl,alfai,'linear','extrap');
    
    %%%%%%%%%%%%%%%%
    % cdi = zeros(size(cli)); % high-Re assumption; flow attached, no skin friction
    cdi = interp1(alfaref,cd,alfai,'linear','extrap');
    %%%%%%%%%%%%%%%%

    % Ensure force coefficients go to zero at ends
    cli(1) = 0; cdi(1) = 0; cli(end) = 0; cdi(end) = 0;
    
    % Calculate planform area
    S = trapz(y,c);  

    % Calculate 
    CL = (1/S) * trapz(y, c.*(cli.*cos(v) - cdi.*sin(v)));
    CD = (1/S) * trapz(y, c.*(cli.*sin(v) + cdi.*cos(v)));

end


%%%%%%%%%%%%%

function [Gammai,vi,alfai,varargout] = calcgamma(yi,ci,thi,AoA,alfa,cl,varargin)

    % Prepare for optional output
    if nargout > 3
        varargout{1} = [];
    end

    % Default settings
    Gamma0 = 0;         % scaling of initial guess for Gamma
    maxiter = 1000;     % maximum iterations
    errtol = 1e-6;      % error tolerance
    relfactor = 0.01;   % relaxation factor    

    % Optional input arguments
    if ~isempty(varargin)
        opts = varargin;
        for ii = 1:length(opts)
            switch opts{ii}
                case 'Gamma0'
                    Gamma0 = opts{ii+1};
                case 'maxiter' 
                    maxiter = opts{ii+1};
                case 'errtol'
                    errtol = opts{ii+1};
                case 'relfactor'
                    relfactor = opts{ii+1};
            end
        end
    end

    % Convert angles to radians
    alfa = pi/180 * alfa;
    AoA = pi/180 * AoA;
    thi = pi/180 * thi;

    % Circulation to be calculated at locations between those in y vector
    yj = 0.5*(yi + circshift(yi,-1));
    yj = yj(1:end-1);

    % Initialize variables
    alfai = zeros(size(yi));
    Gammai = zeros(size(yi));
    cli = zeros(size(yi));

    % Intial guess of Gamma distribution, defined on yi
    Gammai = -Gamma0*(yi-yi(1)).*(yi-yi(end)); % parabolic distribution, zero at ends

    % Ensure end conditions are met
    Gammai(1) = 0;
    Gammai(end) = 0;

    % Iterate until Gamma is converged
    err = 10;
    iter = 0;
    while iter < maxiter && err(end) > errtol
    
        iter = iter + 1;

        % Trailing vortex strengths
        dGammaj = -diff(Gammai);
        
        % Calculate influence matrix
        Yij = 1./( 4* pi * ( yj' - yi ) );
        
        % Calculate induced velocities
        vi = Yij * dGammaj;
        
        % Calculate effective angle of attack
        alfai = AoA + thi - vi;
        %%%%%%%%%%%%%%
        % cli =  interp1(alfa,cl,alfai);
        cli =  interp1(alfa,cl,alfai,'linear','extrap');
        %%%%%%%%%%%%%
        cli(1) = 0; cli(end) = 0;

        % Calculate circulation from lift coefficient
        Gammai_new = cli .* ci/2;
        Gammai_new(1) = 0;
        Gammai_new(end) = 0;

        % Calculate error
        eGamma = Gammai_new - Gammai;
        err(end+1) = max(abs(eGamma));

        % Iterate on Gamma
        Gammai = Gammai + relfactor * eGamma;
    end  

    % Convergence report
    disp(['Total iterations: ',num2str(iter)])

    disp(['err = ', num2str(err(end))])
    if err(end) > errtol
        disp(['NOT CONVERGED WITHIN e = ',num2str(errtol)])
    else
        disp(['CONVERGED WITHIN e = ',num2str(errtol)])
    end
    %%%%%%%%%%
        % Check if angle of attack is out of range
        if any(alfai(2:end-1)>alfa(end)) | any(alfai(2:end-1)<alfa(1))
            disp('ERROR: alfa out of range');
            Gammai = NaN;
            vi = NaN;
            return;
        end
%%%%%%%%%%

    if nargout > 3
        varargout{1} = err;
    end

    % Warning
    % if any(alfai(2:end-1) > alfa(end)) || ...
    %         any(alfai(2:end-1) < alfa(1))
    % 
    %     warning(['Converged effective angle of attack lies outside ', ...
    %          'the supplied airfoil dataset.']);
    % end

    % Convert back to degrees
    alfai = 180/pi * alfai;

end

% function out = checkmonotonic(x)
% 
%     temp = sum(diff(x)<=0);
% 
%     if temp ~= 0
%         disp('ERROR:');
%         out = 0;
%     else
%         out = 1;
%     end
% 
% end

%%%%%%%%%%%%%%%%%
function [flag, indices] = checkmonotonic(x)
%CHECKMONOTONIC Check whether a vector is strictly monotonically increasing.
%
%   [flag, indices] = checkmonotonic(x)
%
%   flag    = 1 if the entire vector is strictly increasing
%             0 otherwise
%
%   indices = indices of the longest contiguous strictly increasing section,
%             with the same row or column orientation as x

    flag = all(diff(x) > 0);

    if flag
        indices = 1:numel(x);
    else
        breaks = find(diff(x) <= 0);

        starts = [1; breaks(:) + 1];
        ends   = [breaks(:); numel(x)];

        [~, iLongest] = max(ends - starts + 1);

        indices = starts(iLongest):ends(iLongest);
    end

    % Match the orientation of x
    if iscolumn(x)
        indices = indices(:);
    end
end
%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%
