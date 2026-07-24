%% test_elliptic_wing.m — Verify LLT solver against classical elliptic wing theory
%
% Tests the liftingline() solver on an untwisted elliptic planform using
% thin-airfoil input data (cl = 2\pi\alpha, cd = 0).
%
% Classical lifting-line theory predicts:
%     dC_L / d\alpha = 2\pi / (1 + 2/AR)                  (lift slope)
%     C_Di = C_L^2 / (\pi AR)                              (induced drag)
%
% Mirrors test_elliptic_wing.py (Python) as the authoritative spec.
%
% Usage:
%     test_elliptic_wing          (run in MATLAB with this directory on path)

clear; clc; close all;

% ═════════════════════════════════════════════════════════════════════════════
%  Parameters
% ═════════════════════════════════════════════════════════════════════════════

AR = 8.0;                           % aspect ratio
N_STATIONS = 31;                    % spanwise discretisation
ALFAS = -10:2:10;                   % angles of attack (deg)
RELAXATION = 0.01;                  % relaxation factor

% Thin-airfoil 2-D lift slope
a0 = 2 * pi;

% Theoretical 3-D lift slope for an elliptic wing (Prandtl)
THEORY_CL_ALPHA_RAD = a0 / (1 + a0 / (pi * AR));
THEORY_CL_ALPHA_DEG = THEORY_CL_ALPHA_RAD * pi / 180;

% Tolerance for passing the test
TOLERANCE_PCT = 2.0;                % lift slope error must be under 2%


% ═════════════════════════════════════════════════════════════════════════════
%  Run solver
% ═════════════════════════════════════════════════════════════════════════════

GEOM_FILE = 'elliptic_AR8.txt';
FORCE_FILE = 'thin_airfoil_data.txt';

fprintf('%s\n', repmat('=', 1, 60));
fprintf('  Lifting-Line Theory --- Elliptic Wing Test\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('\n  Geometry:      elliptic planform, AR = %.0f\n', AR);
fprintf('  Stations:      %d\n', N_STATIONS);
fprintf('  Force data:    thin airfoil (cl = 2\\pi\\alpha, cd = 0)\n');
fprintf('  AoA range:     %.0f to %.0f deg\n', ALFAS(1), ALFAS(end));
fprintf('  Relaxation:    %.2f\n', RELAXATION);
fprintf('\n');

[CL, CD, y, Gamma, v] = liftingline( ...
    GEOM_FILE, FORCE_FILE, ALFAS, 'relfactor', RELAXATION);


% ═════════════════════════════════════════════════════════════════════════════
%  Lift slope check
% ═════════════════════════════════════════════════════════════════════════════

% Linear regression: CL vs alpha (radians)
p = polyfit(deg2rad(ALFAS), CL, 1);
computed_slope_rad = p(1);
computed_intercept = p(2);

error_pct = abs(computed_slope_rad - THEORY_CL_ALPHA_RAD) ...
            / THEORY_CL_ALPHA_RAD * 100;

fprintf('%s\n', repmat('-', 1, 50));
fprintf('  Lift Slope\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('  Computed:      %.4f rad^{-1}  (%.6f deg^{-1})\n', ...
    computed_slope_rad, computed_slope_rad * pi/180);
fprintf('  Theoretical:   %.4f rad^{-1}  (%.6f deg^{-1})\n', ...
    THEORY_CL_ALPHA_RAD, THEORY_CL_ALPHA_DEG);
fprintf('  Error:         %.4f%%\n', error_pct);
fprintf('  Intercept:     %.6f\n', computed_intercept);
if error_pct < TOLERANCE_PCT
    fprintf('  Pass/Fail:     PASS\n');
else
    fprintf('  Pass/Fail:     FAIL\n');
end
fprintf('\n');

% Detailed table
fprintf('%s\n', repmat('-', 1, 50));
fprintf('  Results Table\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%8s %10s %10s %12s\n', 'AoA(deg)', 'CL', 'CD', 'CL^2/(\\pi AR)');
fprintf('%s\n', repmat('-', 1, 42));
for i = 1:length(ALFAS)
    cdp = CL(i)^2 / (pi * AR);
    fprintf('%8d %10.6f %10.6f %12.6f\n', ALFAS(i), CL(i), CD(i), cdp);
end


% ═════════════════════════════════════════════════════════════════════════════
%  Induced drag check
% ═════════════════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('  Induced Drag\n');
fprintf('%s\n', repmat('-', 1, 50));
cd_errors = abs(CD(:) - CL(:).^2 / (pi * AR));
max_cd_error = max(cd_errors);
fprintf('  Max |CD - CL^2/(\\pi AR)| = %.6e\n', max_cd_error);

% Summary
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 60));
if error_pct < TOLERANCE_PCT
    fprintf('  PASSED: Lift slope within %.1f%% of classical theory.\n', TOLERANCE_PCT);
else
    fprintf('  FAILED: Lift slope error %.2f%% exceeds %.1f%%.\n', ...
        error_pct, TOLERANCE_PCT);
end
fprintf('%s\n', repmat('=', 1, 60));