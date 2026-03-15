import tkinter as tk
import pandas as pd
from tkinter import ttk, messagebox, filedialog
from PIL import Image, ImageTk

# Make sure this import matches your backend script name
from sequence_analysis import analyse_sequence, analyse_linear_sequence, DEFAULT_TERMS, DEFAULT_GBLOCK_TAIL


class SequenceApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sequence Analyzer")
        self.geometry("900x600")
        self.resizable(False, False)

        # Standard state variables
        self.internal_exon_data = ""
        self.custom_terms = []
        self.custom_gblock_tail = ""

        # --- NEW VARIABLES FOR ADVANCED SETTINGS ---
        self.pfs_method = tk.StringVar(value="published")
        self.min_dist = tk.IntVar(value=-7)
        self.max_dist = tk.IntVar(value=17)
        self.auto_export = tk.BooleanVar(value=False)

        self.current_bg = None
        self.canvas = tk.Canvas(self, width=900, height=600)
        self.canvas.pack(fill="both", expand=True)

        self.main_menu()

    def set_background(self, image_path):
        try:
            image = Image.open(image_path)
            image = image.resize((900, 600), Image.Resampling.LANCZOS)
            self.bg_image = ImageTk.PhotoImage(image)
            self.canvas.create_image(0, 0, anchor="nw", image=self.bg_image)
        except Exception as e:
            # Fallback if image is missing
            self.canvas.create_rectangle(0, 0, 900, 600, fill="#2E2E2E")

    def clear_canvas(self):
        self.canvas.delete("all")

    def main_menu(self):
        self.clear_canvas()
        self.set_background("images/bg_main.png")

        self.canvas.create_text(450, 80, text="Select Analysis Mode", font=("Helvetica", 26), fill="white")

        self.canvas.create_window(450, 180, window=ttk.Button(self, text="Circular RNA Mode", width=30,
                                                              command=self.circular_mode))
        self.canvas.create_window(450, 240, window=ttk.Button(self, text="Linear Sequence Mode", width=30,
                                                              command=self.linear_mode))
        self.canvas.create_window(450, 300, window=ttk.Button(self, text="Advanced Settings", width=30,
                                                              command=self.advanced_settings))
        self.canvas.create_window(450, 360, window=ttk.Button(self, text="Batch Circular RNA Mode", width=30,
                                                              command=self.batch_circular_mode))
        self.canvas.create_window(450, 420, window=ttk.Button(self, text="Linear Batch Mode", width=30,
                                                              command=self.linear_batch_mode))

    def linear_batch_mode(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV or Excel files", "*.csv *.xlsx")])
        if not file_path:
            return

        try:
            if file_path.endswith(".csv"):
                df = pd.read_csv(file_path)
            elif file_path.endswith(".xlsx"):
                df = pd.read_excel(file_path)
            else:
                messagebox.showerror("Error", "Unsupported file format.")
                return

            if "sequence" not in df.columns:
                messagebox.showerror("Error", "Missing required 'sequence' column.")
                return

            term_list = self.custom_terms if self.custom_terms else DEFAULT_TERMS
            gblock_tail = self.custom_gblock_tail if self.custom_gblock_tail else DEFAULT_GBLOCK_TAIL

            all_results_text = []
            all_dfs = []

            for i, row in df.iterrows():
                name = row.get("name", f"Entry {i + 1}")
                sequence = row["sequence"]

                result_df = analyse_linear_sequence(
                    sequence, term_list, gblock_tail,
                    pfs_method=self.pfs_method.get()
                )

                all_results_text.append(f"=== {name} ===\n{result_df.to_string(index=False)}\n")

                if not result_df.empty:
                    result_df.insert(0, "Sequence_Name", name)
                    all_dfs.append(result_df)

            final_output = "\n\n".join(all_results_text)
            self.show_text_popup("Linear Batch Results", final_output)

            # Auto-export combined dataframe
            if self.auto_export.get() and all_dfs:
                combined_df = pd.concat(all_dfs, ignore_index=True)
                self.export_dataframe(combined_df)

        except Exception as e:
            messagebox.showerror("Error", f"Failed to process file:\n{e}")

    def batch_circular_mode(self):
        self.clear_canvas()
        self.set_background("images/bg_main.png")
        self.canvas.create_text(450, 50, text="Select CSV/Excel file for batch circular RNA analysis",
                                font=("Helvetica", 16), fill="white")

        def open_file():
            file_path = filedialog.askopenfilename(
                title="Select CSV or Excel File",
                filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx *.xls")]
            )
            if not file_path:
                return

            try:
                if file_path.endswith('.csv'):
                    df = pd.read_csv(file_path)
                else:
                    df = pd.read_excel(file_path)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file:\n{e}")
                return

            required_cols = ['name', 'upstream', 'downstream', 'linear_overhang']
            if not all(col in df.columns for col in required_cols):
                messagebox.showerror("Error", f"File must contain columns: {', '.join(required_cols)}")
                return

            term_list = self.custom_terms if self.custom_terms else DEFAULT_TERMS
            gblock_tail = self.custom_gblock_tail if self.custom_gblock_tail else DEFAULT_GBLOCK_TAIL

            final_output = ""
            all_dfs = []

            for index, row in df.iterrows():
                name = str(row['name'])
                upstream = str(row['upstream']).upper()
                downstream = str(row['downstream']).upper()
                linear_overhang = str(row['linear_overhang']).upper()
                internal_exons = str(row['internal_exons']).upper() if 'internal_exons' in df.columns else ""

                circ_results = analyse_sequence(
                    upstream, downstream, term_list, linear_overhang,
                    gblock_tail=gblock_tail,
                    pfs_method=self.pfs_method.get(),
                    min_dist=self.min_dist.get(),
                    max_dist=self.max_dist.get()
                )[-1]  # Only grab the dataframe

                final_output += f"\n=== Circular RNA Results: {name} ===\n"
                final_output += circ_results.to_string(index=False) + "\n"

                if not circ_results.empty:
                    circ_results.insert(0, "Type", "Circular")
                    circ_results.insert(0, "Sequence_Name", name)
                    all_dfs.append(circ_results)

                if internal_exons and internal_exons != 'nan':
                    linear_df = analyse_linear_sequence(
                        internal_exons, term_list, gblock_tail,
                        pfs_method=self.pfs_method.get()
                    )
                    final_output += f"\n--- Internal Exon Linear Analysis: {name} ---\n"
                    final_output += linear_df.to_string(index=False) + "\n"

                    if not linear_df.empty:
                        linear_df.insert(0, "Type", "Internal Exon (Linear)")
                        linear_df.insert(0, "Sequence_Name", name)
                        all_dfs.append(linear_df)

            self.show_text_popup("Batch Circular RNA Results", final_output)

            # Auto-export combined dataframe
            if self.auto_export.get() and all_dfs:
                combined_df = pd.concat(all_dfs, ignore_index=True)
                self.export_dataframe(combined_df)

            self.return_to_menu()

        self.canvas.create_window(450, 100, window=ttk.Button(self, text="Open File", command=open_file))
        self.canvas.create_window(450, 150, window=ttk.Button(self, text="Cancel", command=self.main_menu))

    def circular_mode(self, terms=None):
        self.internal_exon_data = ""
        self.clear_canvas()
        self.set_background("images/bg_upstream.png")

        self.canvas.create_text(450, 50, text="Enter Upstream Sequence", font=("Helvetica", 20), fill="white")
        entry = ttk.Entry(self, width=70)
        self.canvas.create_window(450, 100, window=entry)

        def on_submit_upstream():
            upstream = entry.get().strip().upper()
            if not upstream:
                messagebox.showwarning("Warning", "Please enter the upstream sequence.")
                return
            self.ask_downstream(upstream, terms)

        self.canvas.create_window(400, 150, window=ttk.Button(self, text="Next", command=on_submit_upstream))
        self.canvas.create_window(500, 150, window=ttk.Button(self, text="Cancel", command=self.main_menu))

    def ask_downstream(self, upstream, terms):
        self.clear_canvas()
        self.set_background("images/bg_downstream.png")

        self.canvas.create_text(450, 50, text="Enter Downstream Sequence", font=("Helvetica", 20), fill="white")
        entry = ttk.Entry(self, width=70)
        self.canvas.create_window(450, 100, window=entry)

        def on_submit_downstream():
            downstream = entry.get().strip().upper()
            if not downstream:
                messagebox.showwarning("Warning", "Please enter the downstream sequence.")
                return
            self.ask_linear_overhang(upstream, downstream, terms)

        self.canvas.create_window(400, 150, window=ttk.Button(self, text="Next", command=on_submit_downstream))
        self.canvas.create_window(500, 150, window=ttk.Button(self, text="Cancel", command=self.main_menu))

    def ask_linear_overhang(self, upstream, downstream, terms):
        self.clear_canvas()
        self.set_background("images/bg_linearoverhang.png")

        self.canvas.create_text(450, 50, text="Enter Upstream of Upstream (Linear Overhang)", font=("Helvetica", 18),
                                fill="white")
        entry = ttk.Entry(self, width=70)
        self.canvas.create_window(450, 100, window=entry)

        def on_submit_overhang():
            linearoverhang = entry.get().strip().upper()
            if not linearoverhang:
                messagebox.showwarning("Warning", "Please enter the linear overhang.")
                return
            self.ask_internal_exons(upstream, downstream, linearoverhang, terms)

        self.canvas.create_window(400, 150, window=ttk.Button(self, text="Next", command=on_submit_overhang))
        self.canvas.create_window(500, 150, window=ttk.Button(self, text="Cancel", command=self.main_menu))

    def ask_internal_exons(self, upstream, downstream, linearoverhang, terms):
        self.clear_canvas()
        self.set_background("images/bg_internalexons.png")
        self.canvas.create_text(450, 50, text="Would you like to add internal exons?", font=("Helvetica", 18),
                                fill="white")
        self.canvas.create_window(370, 100, window=ttk.Button(self, text="Yes",
                                                              command=lambda: self.get_internal_exon_loop(upstream,
                                                                                                          downstream,
                                                                                                          linearoverhang,
                                                                                                          terms)))
        self.canvas.create_window(530, 100, window=ttk.Button(self, text="No",
                                                              command=lambda: self.run_circ_analysis(upstream,
                                                                                                     downstream,
                                                                                                     linearoverhang,
                                                                                                     terms)))

    def get_internal_exon_loop(self, upstream, downstream, linearoverhang, terms):
        self.clear_canvas()
        self.set_background("images/bg_internalexons.png")
        self.canvas.create_text(450, 50, text="Enter Internal Exon", font=("Helvetica", 18), fill="white")
        entry = ttk.Entry(self, width=70)
        self.canvas.create_window(450, 100, window=entry)

        def on_submit_exon():
            exon = entry.get().strip().upper()
            if exon:
                self.internal_exon_data += exon
            self.ask_if_last_exon(upstream, downstream, linearoverhang, terms)

        self.canvas.create_window(450, 150, window=ttk.Button(self, text="Submit Exon", command=on_submit_exon))

    def ask_if_last_exon(self, upstream, downstream, linearoverhang, terms):
        self.clear_canvas()
        self.set_background("images/bg_internalexons.png")
        self.canvas.create_text(450, 50, text="Is this the last internal exon?", font=("Helvetica", 18), fill="white")
        self.canvas.create_window(370, 100, window=ttk.Button(self, text="Yes",
                                                              command=lambda: self.run_circ_analysis(upstream,
                                                                                                     downstream,
                                                                                                     linearoverhang,
                                                                                                     terms)))
        self.canvas.create_window(530, 100, window=ttk.Button(self, text="No",
                                                              command=lambda: self.get_internal_exon_loop(upstream,
                                                                                                          downstream,
                                                                                                          linearoverhang,
                                                                                                          terms)))

    def run_circ_analysis(self, upstream, downstream, linearoverhang, terms):
        term_list = terms if terms else (self.custom_terms if self.custom_terms else DEFAULT_TERMS)
        gblock_tail = self.custom_gblock_tail if self.custom_gblock_tail else DEFAULT_GBLOCK_TAIL

        # Run Backend
        circ_results = analyse_sequence(
            upstream, downstream, term_list, linearoverhang,
            gblock_tail=gblock_tail,
            pfs_method=self.pfs_method.get(),
            min_dist=self.min_dist.get(),
            max_dist=self.max_dist.get()
        )

        df = circ_results[-1]
        self.show_results_popup("Circular RNA Results", df)

        # Auto Export Trigger
        if self.auto_export.get() and not df.empty:
            self.export_dataframe(df)

        if self.internal_exon_data:
            linear_df = analyse_linear_sequence(
                self.internal_exon_data, term_list, gblock_tail,
                pfs_method=self.pfs_method.get()
            )
            self.show_results_popup("Internal Exon Linear Analysis", linear_df)

            # Auto Export Trigger for Linear
            if self.auto_export.get() and not linear_df.empty:
                self.export_dataframe(linear_df)

        self.return_to_menu()

    def linear_mode(self, terms=None):
        self.clear_canvas()
        self.set_background("images/bg_linear_mode.png")
        self.canvas.create_text(450, 50, text="Enter Linear Sequence", font=("Helvetica", 20), fill="white")
        entry = ttk.Entry(self, width=80)
        self.canvas.create_window(450, 100, window=entry)

        def on_submit_linear():
            sequence = entry.get().strip().upper()
            if not sequence:
                messagebox.showwarning("Warning", "Please enter a sequence.")
                return

            term_list = terms if terms else (self.custom_terms if self.custom_terms else DEFAULT_TERMS)
            gblock_tail = self.custom_gblock_tail if self.custom_gblock_tail else DEFAULT_GBLOCK_TAIL

            df = analyse_linear_sequence(
                sequence, term_list, gblock_tail,
                pfs_method=self.pfs_method.get()
            )

            self.show_results_popup("Linear Analysis Results", df)

            # Auto Export Trigger
            if self.auto_export.get() and not df.empty:
                self.export_dataframe(df)

            self.return_to_menu()

        self.canvas.create_window(400, 150, window=ttk.Button(self, text="Run Analysis", command=on_submit_linear))
        self.canvas.create_window(500, 150, window=ttk.Button(self, text="Cancel", command=self.main_menu))

    def advanced_settings(self):
        self.clear_canvas()
        self.set_background("images/bg_custom_mode.png")
        self.canvas.create_text(450, 30, text="Advanced Settings", font=("Helvetica", 20), fill="white")

        # Custom Motifs
        self.canvas.create_text(450, 70, text="Custom Motifs (comma-separated):", font=("Helvetica", 12), fill="white")
        motif_entry = ttk.Entry(self, width=80)
        motif_entry.insert(0, ", ".join(self.custom_terms) if self.custom_terms else "")
        self.canvas.create_window(450, 95, window=motif_entry)

        # Custom G-block tail
        self.canvas.create_text(450, 135, text="Custom G-block Tail:", font=("Helvetica", 12), fill="white")
        gblock_entry = ttk.Entry(self, width=80)
        gblock_entry.insert(0, self.custom_gblock_tail if self.custom_gblock_tail else "")
        self.canvas.create_window(450, 160, window=gblock_entry)

        # PFS Scoring Method
        self.canvas.create_text(450, 200, text="PFS Scoring Method:", font=("Helvetica", 12), fill="white")
        pfs_combo = ttk.Combobox(self, textvariable=self.pfs_method, values=["published", "alternative", "average"],
                                 state="readonly", width=30)
        self.canvas.create_window(450, 225, window=pfs_combo)

        # BSJ Distance Filtering
        self.canvas.create_text(450, 265, text="CircRNA Distance Filtering (Min / Max):", font=("Helvetica", 12),
                                fill="white")
        dist_frame = tk.Frame(self, bg="#2E2E2E")  # Dark background to blend in with default empty canvas
        min_entry = ttk.Entry(dist_frame, textvariable=self.min_dist, width=10)
        min_entry.pack(side="left", padx=5)
        max_entry = ttk.Entry(dist_frame, textvariable=self.max_dist, width=10)
        max_entry.pack(side="left", padx=5)
        self.canvas.create_window(450, 290, window=dist_frame)

        # Auto-Export CSV
        export_check = tk.Checkbutton(self, text="Auto-Export Results to CSV", variable=self.auto_export,
                                      font=("Helvetica", 12), bg="white", fg="black")
        self.canvas.create_window(450, 340, window=export_check)

        def on_confirm():
            raw_terms = motif_entry.get().strip().upper()
            self.custom_terms = [t.strip() for t in raw_terms.split(",") if t.strip()]
            self.custom_gblock_tail = gblock_entry.get().strip().upper()

            # Basic validation for distance bounds
            try:
                int(self.min_dist.get())
                int(self.max_dist.get())
                messagebox.showinfo("Settings Saved", "Custom settings applied successfully!")
            except ValueError:
                messagebox.showerror("Error", "Distance boundaries must be integers. Reverting bounds to default.")
                self.min_dist.set(-7)
                self.max_dist.set(17)

        def on_reset():
            self.custom_terms = []
            self.custom_gblock_tail = ""
            motif_entry.delete(0, tk.END)
            gblock_entry.delete(0, tk.END)
            self.pfs_method.set("published")
            self.min_dist.set(-7)
            self.max_dist.set(17)
            self.auto_export.set(False)
            messagebox.showinfo("Settings Reset", "Reverted to default settings.")

        # Action Buttons
        self.canvas.create_window(370, 410, window=ttk.Button(self, text="Confirm", command=on_confirm))
        self.canvas.create_window(530, 410, window=ttk.Button(self, text="Reset to Default", command=on_reset))
        self.canvas.create_window(450, 460, window=ttk.Button(self, text="Back to Menu", command=self.main_menu))

    def custom_mode(self):
        # NOTE: Deprecated mostly by Advanced Settings, but left intact if you still use it directly
        self.clear_canvas()
        self.set_background("images/bg_custom_mode.png")
        self.canvas.create_text(450, 50, text="Enter custom motifs (comma-separated)", font=("Helvetica", 18),
                                fill="white")
        entry = ttk.Entry(self, width=80)
        self.canvas.create_window(450, 100, window=entry)

        def on_submit_custom():
            raw_input = entry.get().strip().upper()
            terms = [t.strip() for t in raw_input.split(",") if t.strip()]
            if not terms:
                messagebox.showwarning("Warning", "Enter at least one motif.")
                return
            self.ask_custom_mode_type(terms)

        self.canvas.create_window(400, 150, window=ttk.Button(self, text="Next", command=on_submit_custom))
        self.canvas.create_window(500, 150, window=ttk.Button(self, text="Cancel", command=self.main_menu))

    def ask_custom_mode_type(self, terms):
        self.clear_canvas()
        self.set_background("images/bg_custom_mode.png")
        self.canvas.create_text(450, 50, text="Run custom terms in which mode?", font=("Helvetica", 18), fill="white")
        self.canvas.create_window(370, 100,
                                  window=ttk.Button(self, text="Circular", command=lambda: self.circular_mode(terms)))
        self.canvas.create_window(530, 100,
                                  window=ttk.Button(self, text="Linear", command=lambda: self.linear_mode(terms)))

    def return_to_menu(self):
        self.internal_exon_data = ""
        self.after(1500, self.main_menu)

    def show_text_popup(self, title, text_content):
        popup = tk.Toplevel(self)
        popup.title(title)
        popup.geometry("800x500")

        text = tk.Text(popup, wrap="none", font=("Courier", 10))
        text.pack(fill="both", expand=True)

        scrollbar = tk.Scrollbar(popup, command=text.yview)
        scrollbar.pack(side="right", fill="y")
        text.config(yscrollcommand=scrollbar.set)

        text.insert("1.0", text_content)
        text.config(state="disabled")

        export_btn = ttk.Button(popup, text="Export Text", command=lambda: self.export_text(text_content))
        export_btn.pack(pady=5, side="bottom")

    def show_results_popup(self, title, df):
        popup = tk.Toplevel(self)
        popup.title(title)
        popup.geometry("800x500")

        frame = tk.Frame(popup)
        frame.pack(fill="both", expand=True, padx=5, pady=5)

        scrollbar = tk.Scrollbar(frame)
        scrollbar.pack(side="right", fill="y")

        text = tk.Text(frame, wrap="none", font=("Courier", 10), yscrollcommand=scrollbar.set)
        text.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=text.yview)

        if df.empty:
            text.insert("1.0", "No matching motifs found.")
        else:
            text.insert("1.0", df.to_string(index=False))

        text.config(state="disabled")

        export_btn = ttk.Button(popup, text="Export Results", command=lambda: self.export_dataframe(df))
        export_btn.pack(pady=5, side="bottom")

    def export_dataframe(self, df):
        if df.empty:
            messagebox.showwarning("No Data", "There is no data to export.")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx")]
        )
        if not file_path:
            return

        try:
            if file_path.endswith(".csv"):
                df.to_csv(file_path, index=False)
            else:
                df.to_excel(file_path, index=False)
            messagebox.showinfo("Export Successful", f"Results saved to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Export Failed", f"Could not save results:\n{e}")

    def export_text(self, text_content):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt")]
        )
        if not file_path:
            return
        try:
            with open(file_path, "w") as f:
                f.write(text_content)
            messagebox.showinfo("Export Successful", f"Results saved to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Export Failed", f"Could not save results:\n{e}")


if __name__ == "__main__":
    app = SequenceApp()
    app.mainloop()