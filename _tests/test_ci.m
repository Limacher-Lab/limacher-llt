%% test_ci.m — CI validation test for the MATLAB lifting-line solver
%
% Runs at 301 spanwise stations on an elliptic AR=8 wing with
% thin-airfoil data and checks against classical theory.
%
% Mirrors _tests/test_ci.py (Python) exactly — same checks, same thresholds.
%
% Exits with code 0 if all checks pass, 1 otherwise.
%
% Usage:
%     test_ci
%     run('_tests/test_ci.m')   (if not in the repo root)

% ── Setup paths ────────────────────────────────────────────────────────────
% Ensure we're in the repo root
test_dir = fileparts(mfilename('fullpath'));  % should be _tests/
repo_dir = fileparts(test_dir);                % repo root
cd(repo_dir);

% ── Parameters ─────────────────────────────────────────────────────────────
AR = 8.0;
STATIONS = 301;
ALFAS = -10:2:10;
FORCE_FILE = fullfile(repo_dir, 'thin_airfoil_data.txt');
RELAXATION = 0.01;

% Theoretical values
a0 = 2 * pi;
theory_slope = a0 / (1 + a0 / (pi * AR));

% ── Generate geometry ──────────────────────────────────────────────────────
y = linspace(-0.5, 0.5, STATIONS)';
c_root = 4.0 / (pi * AR);
c = c_root * sqrt(1.0 - (2.0 * y).^2);
th = zeros(size(y));

geom_file = fullfile(tempdir, 'ci_elliptic.txt');
fid = fopen(geom_file, 'w');
fprintf(fid, 'y,c,th\n');
fprintf(fid, '%.15f,%.15f,%.15f\n', [y, c, th].');
fclose(fid);

% ── Run solver ────────────────────────────────────────────────────────
[CL, CD, ~, ~, ~] = liftingline(geom_file, FORCE_FILE, ALFAS, 'relfactor', RELAXATION);

delete(geom_file);

% ── Checks ─────────────────────────────────────────────────────────────────
exit_code = 0;

% 1. Solver must converge (no NaN)
if any(isnan(CL))
    fprintf('FAIL: Solver returned NaN\n');
    exit_code = 1;
else
    fprintf('PASS: Solver converged (no NaN)\n');
end

% 2. Zero lift at alpha = 0
idx0 = (ALFAS == 0);
cl0 = CL(idx0);
if abs(cl0) > 1e-12
    fprintf('FAIL: CL(alpha=0) = %.3e, expected < 1e-12\n', cl0);
    exit_code = 1;
else
    fprintf('PASS: CL(alpha=0) = %.3e\n', cl0);
end

% 3. Lift slope error
p = polyfit(ALFAS * pi / 180, CL, 1);
slope = p(1);
slope_err = abs(slope - theory_slope) / theory_slope * 100;
if slope_err >= 0.2
    fprintf('FAIL: Lift slope error = %.4f%%, expected < 0.2%%\n', slope_err);
    exit_code = 1;
else
    fprintf('PASS: Lift slope error = %.4f%%\n', slope_err);
end

% 4. Linearity
CL_fit = polyval(p, ALFAS * pi / 180);
residuals = CL(:) - CL_fit(:);
ss_res = sum(residuals.^2);
ss_tot = sum((CL(:) - mean(CL(:))).^2);
r2 = 1 - ss_res / ss_tot;
if r2 <= 0.9999
    fprintf('FAIL: R^2 = %.6f, expected > 0.9999\n', r2);
    exit_code = 1;
else
    fprintf('PASS: R^2 = %.6f\n', r2);
end

% 5. Induced drag
CD_theory = CL.^2 / (pi * AR);
mask = abs(CL) > 1e-10;
if any(mask)
    cd_errors = abs(CD(mask) - CD_theory(mask)) ./ CD_theory(mask) * 100;
    max_cd_err = max(cd_errors);
    if max_cd_err >= 5.0
        fprintf('FAIL: Max induced drag error = %.2f%%, expected < 5%%\n', max_cd_err);
        exit_code = 1;
    else
        fprintf('PASS: Max induced drag error = %.2f%%\n', max_cd_err);
    end
end

% ── Summary ────────────────────────────────────────────────────────────────
fprintf('\n');
if exit_code == 0
    fprintf('All checks passed.\n');
else
    fprintf('Some checks failed.\n');
end

quit(exit_code);