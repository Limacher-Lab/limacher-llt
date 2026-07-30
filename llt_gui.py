"""
llt_gui.py — Limacher Lifting-Line Solver Graphical User Interface

A Windows-friendly tkinter GUI wrapper around the lifting-line solver.
Designed for users who prefer not to edit or run Python scripts directly.

Usage
-----
    python llt_gui.py

Dependencies
------------
    numpy, matplotlib
"""

import sys, io, queue, threading, traceback
from pathlib import Path
import numpy as np

# ══════════════════════════════════════════════════════════════════════════
#  Testable utility functions (no tkinter dependency)
# ══════════════════════════════════════════════════════════════════════════


def parse_geometry_file(path):
    """
    Validate and parse a geometry file.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    y, c, th : ndarray

    Raises
    ------
    FileNotFoundError, ValueError, OSError
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'Geometry file not found: {path}')

    try:
        data = np.loadtxt(path, skiprows=1, delimiter=',')
    except (ValueError, OSError) as e:
        raise ValueError(f'Geometry file could not be read: {e}') from e

    if data.ndim == 1:
        data = data.reshape(-1, 1)  # handle single-column edge case badly→caught below

    if data.shape[1] < 3:
        raise ValueError(
            'Expected at least three comma-separated columns: y, c, th.'
        )

    if data.shape[0] < 2:
        raise ValueError(
            'Geometry file must contain at least 2 spanwise stations.'
        )

    y = data[:, 0]
    c = data[:, 1]
    th = data[:, 2]

    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(c)) or not np.all(np.isfinite(th)):
        raise ValueError(
            'Geometry file contains non-finite values in y, c, or th.'
        )

    if len(y) < 2 or np.any(np.diff(y) <= 0):
        raise ValueError(
            'Spanwise coordinates (y) must be strictly increasing.'
        )

    # Check non-zero span
    if np.abs(y[-1] - y[0]) < 1e-15:
        raise ValueError('Span (y range) must be non-zero.')

    # Negative chord
    if np.any(c < 0):
        raise ValueError('Chord values (c) must not be negative.')

    return y, c, th


def parse_airfoil_file(path):
    """
    Validate and parse an airfoil force-data file.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    alpha, cl, cd : ndarray

    Raises
    ------
    FileNotFoundError, ValueError, OSError
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'Airfoil data file not found: {path}')

    try:
        data = np.loadtxt(path, skiprows=1, delimiter=',')
    except (ValueError, OSError) as e:
        raise ValueError(f'Airfoil data file could not be read: {e}') from e

    if data.ndim == 1:
        data = data.reshape(-1, 1)

    if data.shape[1] < 3:
        raise ValueError(
            'Expected at least three comma-separated columns: Alpha, Cl, Cd.'
        )

    if data.shape[0] < 2:
        raise ValueError(
            'Airfoil data file must contain at least 2 data rows.'
        )

    alpha = data[:, 0]
    cl = data[:, 1]
    cd = data[:, 2]

    if not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(cl)) or not np.all(np.isfinite(cd)):
        raise ValueError(
            'Airfoil data contains non-finite values in Alpha, Cl, or Cd.'
        )

    if np.any(np.diff(alpha) <= 0):
        raise ValueError(
            'Airfoil data angles (Alpha) must be strictly increasing.'
        )

    return alpha, cl, cd


def build_aoa_array(min_aoa, max_aoa, inc):
    """
    Generate an inclusive angle-of-attack array.

    Parameters
    ----------
    min_aoa, max_aoa, inc : float

    Returns
    -------
    ndarray
        Inclusive angle array. Guarantees min_aoa is included.
        max_aoa is included when it falls within a tolerance.
    """
    # Number of steps, including both ends
    n_steps = int(round((max_aoa - min_aoa) / inc)) + 1
    arr = min_aoa + inc * np.arange(n_steps)

    # Clamp to max_aoa within tolerance
    tol = inc * 1e-10
    arr = arr[arr <= max_aoa + tol]

    # Ensure max_aoa is included if close
    if len(arr) > 0 and abs(arr[-1] - max_aoa) > tol and abs(arr[-1] - max_aoa) < inc * 0.5:
        arr = np.append(arr, max_aoa)

    return arr


def validate_angle_sweep(min_aoa, max_aoa, inc):
    """
    Validate angle-of-attack sweep parameters.

    Parameters
    ----------
    min_aoa, max_aoa, inc : float

    Returns
    -------
    str or None
        Error message string, or None if valid.
    """
    if not np.isfinite(min_aoa) or not np.isfinite(max_aoa) or not np.isfinite(inc):
        return 'All angle values must be finite numbers.'

    if inc <= 0:
        return 'Increment must be greater than zero.'

    if max_aoa < min_aoa:
        return 'Maximum angle must be greater than or equal to minimum angle.'

    n = int(round((max_aoa - min_aoa) / inc)) + 1
    if n > 10000:
        return (
            f'The requested sweep contains {n} angles, '
            f'exceeding the 10,000-angle limit.'
        )

    return None  # valid


# ══════════════════════════════════════════════════════════════════════════
#  Queue-based stdout redirector
# ══════════════════════════════════════════════════════════════════════════


class QueueWriter(io.StringIO):
    """
    A file-like object that writes to both an internal StringIO buffer
    and a thread-safe queue for consumption by the GUI main thread.
    """

    def __init__(self, queue):
        super().__init__()
        self._q = queue

    def write(self, s):
        super().write(s)
        self._q.put(s)

    def flush(self):
        pass


# ══════════════════════════════════════════════════════════════════════════
#  GUI application (tkinter)
# ══════════════════════════════════════════════════════════════════════════

# tkinter import deferred so importing llt_gui for tests doesn't require a
# display.  The 'import llt_gui' above only loads the utility functions.
# The GUI class is instantiated only when run via __main__ below.


class LLTGuiApp:
    """
    Main GUI application for the Limacher Lifting-Line Solver.
    """

    def __init__(self, root):
        self.root = root
        self._setup_complete = False

        # Application state
        self.geom_path = None
        self.airfoil_path = None
        self.last_results = None  # (CL, CD, y, Gamma, v, AoA)

        # Thread safety
        self._msg_queue = queue.Queue()
        self._worker_thread = None
        self._sim_running = False
        self._polling = True

        # Defer tkinter imports — they require the event loop
        self._init_ui()

    def _init_ui(self):
        """Set up all tkinter widgets (called after tkinter is imported)."""
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        import matplotlib
        matplotlib.use('TkAgg')

        self.tk = tk
        self.filedialog = filedialog
        self.messagebox = messagebox
        self.FigureCanvasTkAgg = FigureCanvasTkAgg

        self.root.title('Limacher Lifting-Line Solver')
        self.root.minsize(1100, 700)

        # ── Overall layout: 2x2 grid ─────────────────────────────────
        self.root.grid_rowconfigure(0, weight=35)    # top row 35%
        self.root.grid_rowconfigure(1, weight=65)    # bottom row 65%
        self.root.grid_columnconfigure(0, weight=40, uniform='col')   # left 40%
        self.root.grid_columnconfigure(1, weight=60, uniform='col')   # right 60%

        self._build_inputs_panel()
        self._build_geometry_panel()
        self._build_updates_panel()
        self._build_results_panel()

        # Check for default files
        self._load_default_files()

        self._setup_complete = True

    # ── Panel builders ──────────────────────────────────────────────

    def _build_inputs_panel(self):
        """Top-left: Inputs and Run controls."""
        from tkinter import ttk
        frame = ttk.LabelFrame(self.root, text='Inputs', padding=8)
        frame.grid(row=0, column=0, sticky='nsew', padx=4, pady=4)
        frame.grid_rowconfigure(4, weight=1)
        frame.grid_columnconfigure(1, weight=1)

        # Geometry file
        ttk.Label(frame, text='Wing geometry file').grid(
            row=0, column=0, sticky='w', pady=(0, 2))
        self._geom_var = self.tk.StringVar()
        geom_entry = ttk.Entry(frame, textvariable=self._geom_var, state='readonly')
        geom_entry.grid(row=0, column=1, sticky='ew', padx=(2, 4))
        ttk.Button(frame, text='Browse', command=self.browse_geometry_file).grid(
            row=0, column=2, padx=(0, 2))

        # Airfoil data file
        ttk.Label(frame, text='Airfoil data file').grid(
            row=1, column=0, sticky='w', pady=(4, 2))
        self._airfoil_var = self.tk.StringVar()
        af_entry = ttk.Entry(frame, textvariable=self._airfoil_var, state='readonly')
        af_entry.grid(row=1, column=1, sticky='ew', padx=(2, 4))
        ttk.Button(frame, text='Browse', command=self.browse_airfoil_file).grid(
            row=1, column=2, padx=(0, 2))

        # Angle sweep
        sweep_frame = ttk.Frame(frame)
        sweep_frame.grid(row=2, column=0, columnspan=3, sticky='ew', pady=(8, 2))
        ttk.Label(sweep_frame, text='Minimum angle of attack (deg)').grid(
            row=0, column=0, sticky='w', padx=(0, 4))
        self._min_aoa_var = self.tk.StringVar(value='-10')
        ttk.Entry(sweep_frame, textvariable=self._min_aoa_var, width=8).grid(
            row=0, column=1, padx=(0, 12))

        ttk.Label(sweep_frame, text='Maximum angle of attack (deg)').grid(
            row=0, column=2, sticky='w', padx=(0, 4))
        self._max_aoa_var = self.tk.StringVar(value='10')
        ttk.Entry(sweep_frame, textvariable=self._max_aoa_var, width=8).grid(
            row=0, column=3, padx=(0, 12))

        ttk.Label(sweep_frame, text='Increment (deg)').grid(
            row=0, column=4, sticky='w', padx=(0, 4))
        self._inc_var = self.tk.StringVar(value='1')
        ttk.Entry(sweep_frame, textvariable=self._inc_var, width=8).grid(
            row=0, column=5)

        # Angle count display
        self._angle_count_var = self.tk.StringVar(value='Number of requested angles: 21')
        ttk.Label(frame, textvariable=self._angle_count_var).grid(
            row=3, column=0, columnspan=3, sticky='w', pady=(4, 8))

        # Trace angle entries for live count update
        self._min_aoa_var.trace_add('write', self._update_angle_count)
        self._max_aoa_var.trace_add('write', self._update_angle_count)
        self._inc_var.trace_add('write', self._update_angle_count)

        # Run button + status
        run_frame = ttk.Frame(frame)
        run_frame.grid(row=4, column=0, columnspan=3, sticky='ew', pady=(4, 0))
        self._run_btn = ttk.Button(
            run_frame, text='Run simulation', command=self.start_simulation)
        self._run_btn.pack(side='left', padx=(0, 8))
        self._status_var = self.tk.StringVar(value='Ready')
        ttk.Label(run_frame, textvariable=self._status_var).pack(side='left')

    def _build_geometry_panel(self):
        """Bottom-left: Geometry preview."""
        from tkinter import ttk
        from matplotlib.figure import Figure
        frame = ttk.LabelFrame(self.root, text='Geometry preview', padding=4)
        frame.grid(row=1, column=0, sticky='nsew', padx=4, pady=4)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        self._geom_fig = Figure(figsize=(4, 3), dpi=100, tight_layout=True)
        self._geom_ax_chord = self._geom_fig.add_subplot(2, 1, 1)
        self._geom_ax_twist = self._geom_fig.add_subplot(2, 1, 2, sharex=self._geom_ax_chord)

        self._geom_canvas = self.FigureCanvasTkAgg(
            self._geom_fig, master=frame)
        self._geom_canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')

        self._geom_ax_chord.set_ylabel('Chord, c')
        self._geom_ax_twist.set_ylabel('Twist, th (deg)')
        self._geom_ax_twist.set_xlabel('Spanwise coordinate, y')
        self._geom_ax_chord.grid(True, alpha=0.3)
        self._geom_ax_twist.grid(True, alpha=0.3)

    def _build_updates_panel(self):
        """Top-right: Solver updates and errors."""
        from tkinter import ttk
        frame = ttk.LabelFrame(self.root, text='Solver updates', padding=4)
        frame.grid(row=0, column=1, sticky='nsew', padx=4, pady=4)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        self._log_text = self.tk.Text(frame, wrap='word', state='disabled',
                                       height=10, relief='sunken', borderwidth=2)
        self._log_text.grid(row=0, column=0, sticky='nsew')

        scrollbar = ttk.Scrollbar(frame, orient='vertical',
                                   command=self._log_text.yview)
        scrollbar.grid(row=0, column=1, sticky='ns')
        self._log_text['yscrollcommand'] = scrollbar.set

        ttk.Button(frame, text='Clear log', command=self.clear_log).grid(
            row=1, column=0, columnspan=2, pady=(4, 0))

    def _build_results_panel(self):
        """Bottom-right: Final results with notebook tabs."""
        from tkinter import ttk
        from matplotlib.figure import Figure
        frame = ttk.LabelFrame(self.root, text='Results', padding=4)
        frame.grid(row=1, column=1, sticky='nsew', padx=4, pady=4)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        notebook = ttk.Notebook(frame)
        notebook.grid(row=0, column=0, sticky='nsew')

        # Lift curve tab
        lift_frame = ttk.Frame(notebook)
        notebook.add(lift_frame, text='Lift curve')
        lift_frame.grid_rowconfigure(0, weight=1)
        lift_frame.grid_columnconfigure(0, weight=1)

        self._result_fig_lift = Figure(figsize=(5, 3), dpi=100, tight_layout=True)
        self._result_ax_lift = self._result_fig_lift.add_subplot(111)
        self._result_canvas_lift = self.FigureCanvasTkAgg(
            self._result_fig_lift, master=lift_frame)
        self._result_canvas_lift.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        self._result_ax_lift.set_xlabel('Root angle of attack (deg)')
        self._result_ax_lift.set_ylabel('Lift coefficient, CL')
        self._result_ax_lift.grid(True, alpha=0.3)

        # Drag polar tab
        drag_frame = ttk.Frame(notebook)
        notebook.add(drag_frame, text='Drag polar')
        drag_frame.grid_rowconfigure(0, weight=1)
        drag_frame.grid_columnconfigure(0, weight=1)

        self._result_fig_drag = Figure(figsize=(5, 3), dpi=100, tight_layout=True)
        self._result_ax_drag = self._result_fig_drag.add_subplot(111)
        self._result_canvas_drag = self.FigureCanvasTkAgg(
            self._result_fig_drag, master=drag_frame)
        self._result_canvas_drag.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        self._result_ax_drag.set_xlabel('Lift coefficient, CL')
        self._result_ax_drag.set_ylabel('Drag coefficient, CD')
        self._result_ax_drag.grid(True, alpha=0.3)

    # ── Default file loading ────────────────────────────────────────

    def _load_default_files(self):
        """Auto-select example files if found alongside llt_gui.py."""
        script_dir = Path(sys.argv[0] if getattr(sys, 'frozen', False) else __file__).parent
        geom_default = script_dir / 'elliptic_AR8.txt'
        airfoil_default = script_dir / 'thin_airfoil_data.txt'

        if geom_default.exists():
            self._select_geometry_file(str(geom_default))
        if airfoil_default.exists():
            self._select_airfoil_file(str(airfoil_default))

    # ── File browsing ───────────────────────────────────────────────

    def browse_geometry_file(self):
        path = self.filedialog.askopenfilename(
            title='Select wing geometry file',
            filetypes=[('CSV/TXT files', '*.txt *.csv'), ('All files', '*.*')])
        if path:
            self._select_geometry_file(path)

    def browse_airfoil_file(self):
        path = self.filedialog.askopenfilename(
            title='Select airfoil force-data file',
            filetypes=[('CSV/TXT files', '*.txt *.csv'), ('All files', '*.*')])
        if path:
            self._select_airfoil_file(path)

    def _select_geometry_file(self, path):
        self._geom_var.set(path)
        self.geom_path = path
        self.load_and_plot_geometry()

    def _select_airfoil_file(self, path):
        self._airfoil_var.set(path)
        self.airfoil_path = path
        self._validate_and_report_airfoil()

    # ── File validation + geometry plotting ─────────────────────────

    def load_and_plot_geometry(self):
        """Parse geometry file and update geometry preview."""
        path = self._geom_var.get()
        if not path:
            return
        try:
            y, c, th = parse_geometry_file(path)
        except (FileNotFoundError, ValueError, OSError) as e:
            self.append_log(f'Error: {e}')
            # Clear geometry plots
            self._geom_ax_chord.clear()
            self._geom_ax_twist.clear()
            self._update_geom_canvas()
            return

        self.geom_path = path
        self.update_geometry_plots(y, c, th)
        self.append_log(f'Geometry file loaded: {Path(path).name}')

    def _validate_and_report_airfoil(self):
        """Parse airfoil data and log the Alpha range."""
        path = self._airfoil_var.get()
        if not path:
            return
        try:
            alpha, cl, cd = parse_airfoil_file(path)
        except (FileNotFoundError, ValueError, OSError) as e:
            self.append_log(f'Error: {e}')
            return

        self.airfoil_path = path
        self.append_log(
            f'Airfoil data loaded: {Path(path).name}  '
            f'(Alpha range: {alpha[0]:.1f}° to {alpha[-1]:.1f}°, '
            f'{len(alpha)} points)'
        )

    def update_geometry_plots(self, y, c, th):
        """Redraw chord and twist preview plots."""
        self._geom_ax_chord.clear()
        self._geom_ax_twist.clear()

        self._geom_ax_chord.plot(y, c, 'o-', markersize=3, linewidth=1.2)
        self._geom_ax_chord.set_ylabel('Chord, c')
        self._geom_ax_chord.grid(True, alpha=0.3)

        self._geom_ax_twist.plot(y, th, 'o-', markersize=3, linewidth=1.2)
        self._geom_ax_twist.set_ylabel('Twist, th (deg)')
        self._geom_ax_twist.set_xlabel('Spanwise coordinate, y')
        self._geom_ax_twist.grid(True, alpha=0.3)

        self._update_geom_canvas()

    def _update_geom_canvas(self):
        self._geom_canvas.draw_idle()

    # ── Angle count update ──────────────────────────────────────────

    def _update_angle_count(self, *_):
        """Update the 'Number of requested angles' label."""
        try:
            min_aoa = float(self._min_aoa_var.get())
            max_aoa = float(self._max_aoa_var.get())
            inc = float(self._inc_var.get())
            if inc > 0 and max_aoa >= min_aoa:
                arr = build_aoa_array(min_aoa, max_aoa, inc)
                self._angle_count_var.set(f'Number of requested angles: {len(arr)}')
                return
        except (ValueError, ZeroDivisionError):
            pass
        self._angle_count_var.set('Number of requested angles: —')

    # ── Logging ─────────────────────────────────────────────────────

    def append_log(self, message):
        """Append a message to the solver updates text widget."""
        self._log_text.configure(state='normal')
        self._log_text.insert('end', message + '\n')
        self._log_text.see('end')
        self._log_text.configure(state='disabled')
        self.root.update_idletasks()

    def clear_log(self):
        self._log_text.configure(state='normal')
        self._log_text.delete('1.0', 'end')
        self._log_text.configure(state='disabled')

    # ── Simulation ──────────────────────────────────────────────────

    def start_simulation(self):
        """Validate inputs and start the solver in a worker thread."""
        if self._sim_running:
            return

        # Validate geometry file
        geom_path = self._geom_var.get()
        try:
            parse_geometry_file(geom_path)
        except (FileNotFoundError, ValueError, OSError) as e:
            self.append_log(f'Error: {e}')
            return

        # Validate airfoil file
        airfoil_path = self._airfoil_var.get()
        try:
            parse_airfoil_file(airfoil_path)
        except (FileNotFoundError, ValueError, OSError) as e:
            self.append_log(f'Error: {e}')
            return

        # Validate angle sweep
        try:
            min_aoa = float(self._min_aoa_var.get())
            max_aoa = float(self._max_aoa_var.get())
            inc = float(self._inc_var.get())
        except ValueError:
            self.append_log('Error: All angle values must be valid numbers.')
            return

        err = validate_angle_sweep(min_aoa, max_aoa, inc)
        if err is not None:
            self.append_log(f'Error: {err}')
            return

        # Everything valid — build AoA array
        AoA = build_aoa_array(min_aoa, max_aoa, inc)

        # Clear previous result plots
        self._result_ax_lift.clear()
        self._result_ax_lift.set_xlabel('Root angle of attack (deg)')
        self._result_ax_lift.set_ylabel('Lift coefficient, CL')
        self._result_ax_lift.grid(True, alpha=0.3)
        self._result_canvas_lift.draw_idle()

        self._result_ax_drag.clear()
        self._result_ax_drag.set_xlabel('Lift coefficient, CL')
        self._result_ax_drag.set_ylabel('Drag coefficient, CD')
        self._result_ax_drag.grid(True, alpha=0.3)
        self._result_canvas_drag.draw_idle()

        # Log start
        self.append_log(
            f'\n--- Simulation: AoA ({min_aoa}° to {max_aoa}°, '
            f'inc {inc}°, {len(AoA)} angles) ---'
        )

        # UI state
        self._sim_running = True
        self._run_btn.configure(state='disabled')
        self._status_var.set('Running')
        self._polling = True

        # Set up message queue and stdout redirect
        self._msg_queue = queue.Queue()
        self._captured_stdout = io.StringIO()

        # Run solver in thread
        self._worker_thread = threading.Thread(
            target=self._simulation_worker,
            args=(geom_path, airfoil_path, AoA),
            daemon=True
        )
        self._worker_thread.start()

        # Start polling the message queue
        self._poll_queue()

    def _simulation_worker(self, geom_path, airfoil_path, AoA):
        """Run the solver (worker thread)."""
        from liftingline import liftingline

        # Redirect stdout to queue
        old_stdout = sys.stdout
        sys.stdout = QueueWriter(self._msg_queue)

        try:
            CL, CD, y, Gamma, v = liftingline(
                geom_path, airfoil_path, AoA=AoA, relfactor=0.01
            )
            # Put result back via queue
            self._msg_queue.put(('__RESULT__', (CL, CD, y, Gamma, v, AoA)))
        except Exception as exc:
            self._msg_queue.put(('__ERROR__', (type(exc).__name__, str(exc),
                                               traceback.format_exc())))
        finally:
            sys.stdout = old_stdout
            self._msg_queue.put('__DONE__')

    def _poll_queue(self):
        """Poll the message queue from the main thread (called via root.after)."""
        if not self._polling:
            return

        try:
            while True:
                msg = self._msg_queue.get_nowait()

                if isinstance(msg, tuple):
                    if msg[0] == '__RESULT__':
                        self._handle_simulation_success(*msg[1])
                        continue
                    elif msg[0] == '__ERROR__':
                        self._handle_simulation_failure(*msg[1])
                        continue
                elif msg == '__DONE__':
                    continue

                # Regular stdout text
                self.append_log(msg.rstrip())

        except queue.Empty:
            pass

        # Check if thread is still alive
        if self._worker_thread is not None and self._worker_thread.is_alive():
            self.root.after(100, self._poll_queue)
        else:
            self._finalize_run()

    def _finalize_run(self):
        """Clean up after simulation thread finishes."""
        self._sim_running = False
        self._worker_thread = None
        self._run_btn.configure(state='normal')
        self.root.update_idletasks()

    def _handle_simulation_success(self, CL, CD, y, Gamma, v, AoA):
        """Process successful simulation results."""
        self.last_results = (CL, CD, y, Gamma, v, AoA)

        # Check for invalid (NaN) points
        valid_mask = ~(np.isnan(CL) | np.isnan(CD))
        n_nan = np.sum(~valid_mask)

        if n_nan > 0:
            nan_angles = AoA[~valid_mask]
            self.append_log(
                f'\nWarning: {n_nan} invalid point(s) at AoA = '
                f'{", ".join(f"{a:.2f}°" for a in nan_angles)}'
            )
            self._status_var.set('Completed with invalid points')
        else:
            self._status_var.set('Completed')
            self.append_log('\nSimulation completed successfully.')

        # Plot results
        self._result_ax_lift.clear()
        self._result_ax_lift.plot(AoA, CL, 'o-', markersize=4, linewidth=1.2)
        self._result_ax_lift.set_xlabel('Root angle of attack (deg)')
        self._result_ax_lift.set_ylabel('Lift coefficient, CL')
        self._result_ax_lift.grid(True, alpha=0.3)
        self._result_ax_lift.autoscale()
        self._result_canvas_lift.draw_idle()

        self._result_ax_drag.clear()
        self._result_ax_drag.plot(CL, CD, 'o-', markersize=4, linewidth=1.2)
        self._result_ax_drag.set_xlabel('Lift coefficient, CL')
        self._result_ax_drag.set_ylabel('Drag coefficient, CD')
        self._result_ax_drag.grid(True, alpha=0.3)
        self._result_ax_drag.autoscale()
        self._result_canvas_drag.draw_idle()

    def _handle_simulation_failure(self, exc_type, exc_msg, exc_tb):
        """Process simulation failure."""
        self.append_log(f'\nSimulation failed: {exc_type} — {exc_msg}')
        self.append_log(f'Traceback:\n{exc_tb}')
        self._status_var.set('Failed')

        # Clear result axes
        self._result_ax_lift.clear()
        self._result_ax_lift.set_xlabel('Root angle of attack (deg)')
        self._result_ax_lift.set_ylabel('Lift coefficient, CL')
        self._result_ax_lift.grid(True, alpha=0.3)
        self._result_canvas_lift.draw_idle()

        self._result_ax_drag.clear()
        self._result_ax_drag.set_xlabel('Lift coefficient, CL')
        self._result_ax_drag.set_ylabel('Drag coefficient, CD')
        self._result_ax_drag.grid(True, alpha=0.3)
        self._result_canvas_drag.draw_idle()

    def on_close(self):
        """Clean shutdown."""
        self._polling = False
        self.root.destroy()


# ══════════════════════════════════════════════════════════════════════════
#  Entry point
# ══════════════════════════════════════════════════════════════════════════

def main():
    import tkinter as tk
    root = tk.Tk()
    app = LLTGuiApp(root)
    root.protocol('WM_DELETE_WINDOW', app.on_close)
    root.mainloop()


if __name__ == '__main__':
    main()