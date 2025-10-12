import tkinter as tk
from tkinter import ttk, StringVar
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math

class DataPlotterApp:
    """
    A Tkinter application featuring self.number_of_params interactive data control parameters
    and a real-time Matplotlib plot.
    Activation of an parameter is handled via a mouse click on the parameter.
    """
    def __init__(self, root):
        self.root = root
        self.root.title("Interactive Data Plotter")
        self.root.state('zoomed')  # Maximize the window

        # Internal state for the self.number_of_params parameters: value, delta, and Tkinter variables
        self.parameters = []
        self.parameter_frames = []
        self.active_parameter_index = 0  # First parameter is active by default
        self.param_rows = 2
        self.params_per_row = 7
        self.number_of_params = self.param_rows * self.params_per_row

        self._initialize_grid_weights()
        self._initialize_parameter_data()
        self._create_widgets()
        self._setup_bindings()

        self.update_plot()  # Initial plot
        self.activate_parameter(self.active_parameter_index)  # Initial highlighting

    def _initialize_grid_weights(self):
        """Set grid weights for main window to achieve proportional sizing."""
        # Main layout: 3 rows (8%, 8%, 84%) and self.params_per_row columns
        total_weight = 100
        parameter_row_weight = math.ceil(total_weight * 0.08)  # Calculate a suitable weight for the parameter rows (8% of total height)
        self.root.grid_rowconfigure(0, weight=parameter_row_weight)  # Row 0 (Parameters 1-self.params_per_row): 8%
        self.root.grid_rowconfigure(1, weight=parameter_row_weight)  # Row 1 (Parameters 7-self.number_of_params): 8%
        self.root.grid_rowconfigure(2, weight=total_weight - 2 * parameter_row_weight)  # Row 2 (Plot): 84%
        for i in range(self.params_per_row):
            self.root.grid_columnconfigure(i, weight=1)  # self.params_per_row columns for the self.params_per_row parameters per row

    def _initialize_parameter_data(self):
        """Initialize internal data model and Tkinter variables for all self.number_of_params parameters."""
        for i in range(self.number_of_params):
            parameter_data = {
                'value': 100.0,
                'delta': 1.0,
                'value_var': StringVar(value=f"{100.0:.6e}"),
                'delta_var': StringVar(value=f"{1.0:.6e}")
            }
            self.parameters.append(parameter_data)

    def _create_widgets(self):
        """
        Create and place all self.number_of_params control parameters and the Matplotlib plot.
        The activation is triggered by a mouse click on the parameter frame or its contents.
        """
        # 1. Control Parameters (2 rows x self.params_per_row columns)
        for i in range(self.number_of_params):
            parameter_num = i + 1
            parameter = self.parameters[i]
            frame = ttk.Frame(self.root, padding="5 5 5 5", relief=tk.RIDGE)  # Create a Frame for the parameter (for grouping and highlighting)
            row = i // self.params_per_row  # Calculate grid position: Row 0/1, Column 0-5
            col = i % self.params_per_row

            frame.grid(row=row, column=col, sticky="nsew", padx=5, pady=5)
            self.parameter_frames.append(frame)
            frame.bind('<Button-1>', lambda event, idx=i: self.activate_parameter(idx))  # Bind mouse click event to the frame to activate the parameter

            for j in range(3):  # Configure frame's internal grid (3 rows: Title, Value, Delta)
                frame.grid_rowconfigure(j, weight=1)
            frame.grid_columnconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=2)
            # A. Title Label (Row 0)
            title_label = ttk.Label(frame, text=f"Parameter {parameter_num}", font=('Inter', 10, 'bold'))
            title_label.grid(row=0, column=0, columnspan=2, sticky='w', pady=(0, 5))
            # B. Value Entry (Row 1)
            value_label = ttk.Label(frame, text="Value:", font=('Inter', 9))
            value_label.grid(row=1, column=0, sticky='w')
            value_entry = ttk.Entry(frame, textvariable=parameter['value_var'], width=10, state='readonly')
            value_entry.grid(row=1, column=1, sticky='ew', padx=(5, 0))
            # C. Delta Entry (Row 2)
            delta_label = ttk.Label(frame, text="Delta:", font=('Inter', 9))
            delta_label.grid(row=2, column=0, sticky='w')
            delta_entry = ttk.Entry(frame, textvariable=parameter['delta_var'], width=10, state='readonly')
            delta_entry.grid(row=2, column=1, sticky='ew', padx=(5, 0))
            for widget in (title_label, value_label, value_entry, delta_label, delta_entry):  # Bind child widgets to the same activation command so clicking on them works too
                widget.bind('<Button-1>', lambda event, idx=i: self.activate_parameter(idx))

        # Matplotlib Plot (Row 2)
        plot_frame = ttk.Frame(self.root, padding="10 10 10 10", relief=tk.SUNKEN)
        plot_frame.grid(row=2, column=0, columnspan=self.params_per_row, sticky="nsew", padx=5, pady=5)
        plot_frame.grid_rowconfigure(0, weight=1)
        plot_frame.grid_columnconfigure(0, weight=1)

        # Create Matplotlib Figure
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Parameter 1 (X) vs. Parameter 2 (Y) Plot")
        self.ax.set_xlabel("Parameter 1 Value")
        self.ax.set_ylabel("Parameter 2 Value")
        self.scatter, = self.ax.plot([], [], 'o', color='royalblue', markersize=10)
        self.ax.grid(True, linestyle='--', alpha=0.6)

        # Embed the figure into Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.grid(row=0, column=0, sticky="nsew")

    def _setup_bindings(self):
        """Bind arrow keys to the root window."""
        # Ensure that the arrow keys only work when the main window is focused.
        self.root.bind('<Up>', self.handle_key)
        self.root.bind('<Down>', self.handle_key)
        self.root.bind('<Left>', self.handle_key)
        self.root.bind('<Right>', self.handle_key)

    def activate_parameter(self, index):
        """Sets the currently active parameter, visually highlights it, and ensures key focus."""
        if self.active_parameter_index is not None and self.active_parameter_index < len(self.parameter_frames):
            self.parameter_frames[self.active_parameter_index]['style'] = 'TFrame'  # Reset color of the previously active frame
        self.active_parameter_index = index
        self.root.style = ttk.Style()  # Define a style for the active frame and apply it
        self.root.style.configure('Active.TFrame', background='#e0ffe0', relief='solid', borderwidth=2)  # Use a solid border with a light green background for clear highlighting
        self.parameter_frames[self.active_parameter_index]['style'] = 'Active.TFrame'
        self.root.focus_set()  # Crucially, set focus back to the root window so arrow key bindings work
        print(f"Parameter {index + 1} activated. Use arrow keys to adjust its value.")

    def handle_key(self, event):
        """Processes arrow key presses to modify the active parameter's value or delta."""
        if self.active_parameter_index is None:
            return
        parameter = self.parameters[self.active_parameter_index]
        key = event.keysym
        if key == 'Up':
            parameter['delta'] *= 2
        elif key == 'Down':
            parameter['delta'] /= 2
        elif key == 'Left':
            parameter['value'] -= parameter['delta']
        elif key == 'Right':
            parameter['value'] += parameter['delta']

        # Update the UI and the plot after every change
        self.update_entry_fields()
        self.update_plot()

    def update_entry_fields(self):
        """Update the Tkinter entry variables from the internal data model."""
        for parameter in self.parameters:
            # Update the StringVars with formatted (2 decimal places) float values
            parameter['value_var'].set(f"{parameter['value']:.6e}")
            parameter['delta_var'].set(f"{parameter['delta']:.6e}")

    def update_plot(self):
        """Update the Matplotlib plot with the new coordinates (Parameter1.value, Parameter2.value)."""
        x_val = self.parameters[0]['value']  # Parameter 1 value
        y_val = self.parameters[1]['value']  # Parameter 2 value

        self.scatter.set_data([x_val], [y_val])  # Update scatter plot data

        # Adjust plot limits dynamically (10% padding based on current values)
        # Consider a base range to prevent initial plot from being too small
        base_range = 10
        all_values = [a['value'] for a in self.parameters]

        # Calculate overall min/max value across all self.number_of_params parameters
        min_v = min(all_values)
        max_v = max(all_values)

        # Determine plot limits based on the overall range plus padding
        padding = (max_v - min_v) * 0.1
        if padding < base_range:
            padding = base_range

        x_min = min_v - padding
        x_max = max_v + padding
        y_min = min_v - padding
        y_max = max_v + padding

        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

        # Redraw the canvas
        self.canvas.draw_idle()

class CurveSimGUIfit:
    def __init__(self, p, bodies, time_s0, time_d, measured_tt):
        root = tk.Tk()
        app = DataPlotterApp(root)
        root.focus_set()  # Give focus to the main window so key presses are immediately captured
        root.mainloop()

    def save_lmfit_results(self, p):
        pass
