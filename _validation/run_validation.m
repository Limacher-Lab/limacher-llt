%% run_validation.m — Run full LLT validation and save results to CSV
%
% Does three things:
% 1. Convergence study: lift slope error vs number of spanwise stations
%    -> Output: validation_matlab_convergence.csv
% 2. Lift curve and drag polar at 301 stations
%    -> Output: validation_matlab_data.csv
% 3. Prints a summary to console
%
% Mirrors _validation/run_validation.py (Python).
% Compatible with both MATLAB and Octave.
%
% Usage:
%     run_validation            (from repo root)
%     run('_validation/run_validation.m')

% ── Setup paths ────────────────────────────────────────────────────────────
test_dir = fileparts(mfilename('fullpath'));  % _validation/
repo_dir = fileparts(test_dir);                % repo root
cd(repo_dir);

% ── Parameters ─────────────────────────────────────────────────────────────
AR = 8.0;
ALFAS = -10:2:10;
RELAXATION = 0.01;
STATIONS = [31, 151, 301, 601];        % for convergence study
N_DETAIL = 301;                                % for detailed lift/drag data

a0 = 2 * pi;
theory_slope = a0 / (1 + a0 / (pi * AR));
force_file = fullfile(repo_dir, 'thin_airfoil_data.txt');


% ═════════════════════════════════════════════════════════════════════════════
%  1. Convergence study
% ═════════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('=', 1, 60));
fprintf('  LLT Validation --- MATLAB Implementation\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('\n  AR = %.0f, AoA = %.0f to %.0f deg, relfactor = %.2f\n', ...
    AR, ALFAS(1), ALFAS(end), RELAXATION);
fprintf('\n--- Convergence Study ---\n');
fprintf('%10s  %16s  %10s  %10s\n', 'Stations', 'dCL/da (rad-1)', 'Error (%)', 'Avg iters');
fprintf('%s\n', repmat('-', 1, 52));

conv_data = [];
for s = 1:length(STATIONS)
    n = STATIONS(s);

    % Generate elliptic geometry
    geom_file = fullfile(tempdir, sprintf('elliptic_%d.txt', n));
    create_elliptic_wing(geom_file, n, AR);

    % Run solver
    [CL, CD, ~, ~, ~] = liftingline(geom_file, force_file, ALFAS, 'relfactor', RELAXATION);

    delete(geom_file);

    % Compute lift slope and error
    if any(isnan(CL))
        slope = NaN;
        err = NaN;
        fprintf('%10d  %16s  %10s  %10.1f  (did not converge)\n', n, '-', '-', NaN);
    else
        p = polyfit(ALFAS * pi / 180, CL, 1);
        slope = p(1);
        err = abs(slope - theory_slope) / theory_slope * 100;
        fprintf('%10d  %16.6f  %9.4f%%  %10.1f\n', n, slope, err, NaN);
    end

    conv_data = [conv_data; n, slope, err, NaN];
end

% Save convergence CSV
conv_csv = fullfile(test_dir, 'validation_matlab_convergence.csv');
fid = fopen(conv_csv, 'w');
fprintf(fid, 'Stations,dCL_dalpha,Error_pct,Avg_iters\n');
for i = 1:size(conv_data, 1)
    fprintf(fid, '%d,%.6f,%.4f,%.1f\n', ...
        conv_data(i,1), conv_data(i,2), conv_data(i,3), conv_data(i,4));
end
fclose(fid);
fprintf('\n  Saved: %s\n', conv_csv);


% ═════════════════════════════════════════════════════════════════════════════
%  2. Detailed lift-curve and drag-polar data (at N_DETAIL stations)
% ═════════════════════════════════════════════════════════════════════════════

fprintf('\n--- Detailed Results (%d stations) ---\n', N_DETAIL);

geom_file = fullfile(tempdir, 'elliptic_detail.txt');
create_elliptic_wing(geom_file, N_DETAIL, AR);

[CL, CD, ~, ~, ~] = liftingline(geom_file, force_file, ALFAS, 'relfactor', RELAXATION);
delete(geom_file);

CL_theory = theory_slope * ALFAS * pi / 180;
CD_theory = CL_theory.^2 / (pi * AR);

fprintf('%8s %10s %10s %10s %10s\n', 'AoA(deg)', 'CL', 'CD', 'CL_theory', 'CD_theory');
fprintf('%s\n', repmat('-', 1, 50));
for i = 1:length(ALFAS)
    fprintf('%8d %10.6f %10.6f %10.6f %10.6f\n', ...
        ALFAS(i), CL(i), CD(i), CL_theory(i), CD_theory(i));
end

% Save detailed data CSV
data_csv = fullfile(test_dir, 'validation_matlab_data.csv');
fid = fopen(data_csv, 'w');
fprintf(fid, 'AoA,CL,CD,CL_theory,CD_theory\n');
for i = 1:length(ALFAS)
    fprintf(fid, '%d,%.6f,%.6f,%.6f,%.6f\n', ...
        ALFAS(i), CL(i), CD(i), CL_theory(i), CD_theory(i));
end
fclose(fid);
fprintf('\n  Saved: %s\n', data_csv);

% ── Summary ─────────────────────────────────────────────────────────────────
fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('  Validation complete.\n');
fprintf('%s\n', repmat('=', 1, 60));