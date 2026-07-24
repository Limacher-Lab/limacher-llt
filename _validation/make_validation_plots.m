%% make_validation_plots.m — Generate validation plots for the LLT solvers
%
% Runs the MATLAB lifting-line solver on an elliptic AR=8 wing and
% produces three figures: lift curve, drag polar, and convergence study.
%
% Saves PNGs to the _validation/ directory, overwriting existing files.
% Compatible with both MATLAB and Octave.
%
% Mirrors _validation/make_validation_plots.py (Python).
%
% Usage:
%     make_validation_plots
%     run('_validation/make_validation_plots.m')

% ── Setup paths ────────────────────────────────────────────────────────────
test_dir = fileparts(mfilename('fullpath'));  % _validation/
repo_dir = fileparts(test_dir);                % repo root
cd(repo_dir);

% ── Parameters ─────────────────────────────────────────────────────────────
AR = 8.0;
ALFAS = -10:2:10;

% Run solver
[CL, CD, ~, ~, ~] = liftingline( ...
    'elliptic_AR8.txt', 'thin_airfoil_data.txt', ...
    ALFAS, 'relfactor', 0.01);

% Theoretical values
a0 = 2 * pi;
theory_slope = a0 / (1 + a0 / (pi * AR));
theory_CL = theory_slope * ALFAS * pi / 180;
theory_CDi = theory_CL.^2 / (pi * AR);

% Compute lift slope error
p = polyfit(ALFAS * pi / 180, CL, 1);
slope = p(1);
err = abs(slope - theory_slope) / theory_slope * 100;

% Set font size defaults
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontSize', 12);

% ═════════════════════════════════════════════════════════════════════════════
%  Figure 1: Lift curve
% ═════════════════════════════════════════════════════════════════════════════
figure('Visible', 'off', 'Position', [100, 100, 560, 420]);
plot(ALFAS, CL, 'bo-', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(ALFAS, theory_CL, 'r--', 'LineWidth', 2);
xlabel('Angle of Attack (deg)');
ylabel('C_L');
title(sprintf('Elliptic Wing (AR=%.0f) - Lift Curve', AR));
legend({'LLT (computed)', sprintf('Theory: dC_L/d\\alpha = %.4f rad^{-1}', theory_slope)}, ...
    'Location', 'northwest');
grid on;
% Add error text box
text(0.97, 0.05, sprintf('Error: %.3f%%', err), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'BackgroundColor', [1, 0.96, 0.8], 'EdgeColor', 'black', ...
    'FontSize', 12);
saveas(gcf, fullfile(test_dir, 'validation_lift_curve.png'));
close;

% ═════════════════════════════════════════════════════════════════════════════
%  Figure 2: Drag polar
% ═════════════════════════════════════════════════════════════════════════════
figure('Visible', 'off', 'Position', [100, 100, 560, 420]);
plot(CL, CD, 'bo-', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(theory_CL, theory_CDi, 'r--', 'LineWidth', 2);
xlabel('C_L');
ylabel('C_D');
title(sprintf('Elliptic Wing (AR=%.0f) - Drag Polar', AR));
legend({'LLT (computed)', 'Theory: C_{Di} = C_L^2 / (\\pi AR)'}, ...
    'Location', 'northwest');
grid on;
saveas(gcf, fullfile(test_dir, 'validation_drag_polar.png'));
close;

% ═════════════════════════════════════════════════════════════════════════════
%  Figure 3: Convergence study
% ═════════════════════════════════════════════════════════════════════════════
stations = [31, 101, 301, 601];
errors   = [0.6143, 0.1689, 0.1073, 0.0756];

figure('Visible', 'off', 'Position', [100, 100, 560, 420]);
semilogx(stations, errors, 'go-', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
plot(xlim, [0.1, 0.1], 'r--', 'LineWidth', 1.5);
xlabel('Number of Spanwise Stations');
ylabel('Lift Slope Error (%)');
title('Convergence Study - Error vs Discretization');
legend({'LLT (computed)', '0.1% threshold'}, 'Location', 'northeast');
grid on;
set(gca, 'XTick', stations);
set(gca, 'XTickLabel', arrayfun(@num2str, stations, 'UniformOutput', false));
text(0.97, 0.95, '601 stations: 0.076%', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'BackgroundColor', [0.9, 1, 0.9], 'EdgeColor', 'black', ...
    'FontSize', 12);
saveas(gcf, fullfile(test_dir, 'validation_convergence.png'));
close;

fprintf('Done\n');