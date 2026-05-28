#!/usr/bin/env python3
"""
CircHIT GUI
-----------
Design gRNA target sites and scan for genome‑wide off‑targets.
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import threading
import os
import pandas as pd
from PIL import Image, ImageTk

# Import backend functions
from circ_hit_backend import (
    analyse_sequence,
    analyse_linear_sequence,
    load_fd_scores,
    get_24nt_guide,
    run_off_target_for_guide,
    DEFAULT_TERMS,
    DEFAULT_GBLOCK_TAIL,
)


class SequenceApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("CircHIT – Circular & Linear gRNA Designer with Off‑Target Scanner")
        self.geometry("900x650")
        self.resizable(False, False)

        # Internal state
        self.internal_exon_data = ""
        self.custom_terms = []
        self.custom_gblock_tail = ""
        self.fd_csv_path = tk.StringVar(value="00_inputs/Raw_PFS_FD.csv")
        self.min_dist = tk.IntVar(value=-7)
        self.max_dist = tk.IntVar(value=17)
        self.auto_export = tk.BooleanVar(value=False)
        self.top_n = tk.IntVar(value=1000)
        self.mm_cutoff = tk.IntVar(value=5)

        # FD scores loaded once
        self.fd_map = {}
        self.load_fd_scores()

        self.current_bg = None
        self.canvas = tk.Canvas(self, width=900, height=650)
        self.canvas.pack(fill="both", expand=True)

        self.main_menu()

    # ---------- image helper ----------
    def set_background(self, image_path):
        try:
            image = Image.open(image_path)
            image = image.resize((900, 650), Image.Resampling.LANCZOS)
            self.bg_image = ImageTk.PhotoImage(image)
            self.canvas.create_image(0, 0, anchor="nw", image=self.bg_image)
        except Exception:
            self.canvas.create_rectangle(0, 0, 900, 650, fill="#2E2E2E")

    def clear_canvas(self):
        self.canvas.delete("all")

    # ---------- FD loading ----------
    def load_fd_scores(self):
        csv_path = self.fd_csv_path.get()
        if os.path.exists(csv_path):
            self.fd_map = load_fd_scores(csv_path)
        else:
            self.fd_map = {}
            print("[WARNING] FD CSV not found; scores will be 0.")

    # ---------- Menu ----------
    def main_menu(self):
        self.clear_canvas()
        self.set_background("images/bg_main.png")
        self.canvas.create_text(450, 60, text="CircHIT", font=("Helvetica", 32, "bold"), fill="white")
        self.canvas.create_text(450, 110, text="gRNA Design & Off‑Target Scanner", font=("Helvetica", 14), fill="lightgray")

        y_start = 180
        dy = 55
        self.canvas.create_window(450, y_start, window=ttk.Button(self, text="Circular RNA Mode", width=30,
                                                                  command=self.circular_mode))
        self.canvas.create_window(450, y_start+dy, window=ttk.Button(self, text="Linear Sequence Mode", width=30,
                                                                     command=self.linear_mode))
        self.canvas.create_window(450, y_start+2*dy, window=ttk.Button(self, text="Off‑Target Scan", width=30,
                                                                       command=self.off_target_window))
        self.canvas.create_window(450, y_start+3*dy, window=ttk.Button(self, text="Batch Circular RNA Mode", width=30,
                                                                       command=self.batch_circular_mode))
        self.canvas.create_window(450, y_start+4*dy, window=ttk.Button(self, text="Linear Batch Mode", width=30,
                                                                       command=self.linear_batch_mode))
        self.canvas.create_window(450, y_start+5*dy, window=ttk.Button(self, text="Advanced Settings", width=30,
                                                                       command=self.advanced_settings))

    # ---------- Advanced Settings ----------
    def advanced_settings(self):
        self.clear_canvas()
        self.set_background("images/bg_custom_mode.png")
        self.canvas.create_text(450, 30, text="Advanced Settings", font=("Helvetica", 20), fill="white")

        # FD CSV path
        self.canvas.create_text(450, 70, text="FD Score CSV:", font=("Helvetica", 12), fill="white")
        fd_entry = ttk.Entry(self, textvariable=self.fd_csv_path, width=70)
        self.canvas.create_window(450, 95, window=fd_entry)
        ttk.Button(self, text="Browse", command=lambda: self.browse_fd_csv(fd_entry)).place(x=760, y=85)

        # Custom motifs
        self.canvas.create_text(450, 135, text="Custom Motifs (comma-separated):", font=("Helvetica", 12), fill="white")
        motif_entry = ttk.Entry(self, width=80)
        motif_entry.insert(0, ", ".join(self.custom_terms) if self.custom_terms else "")
        self.canvas.create_window(450, 160, window=motif_entry)

        # Custom g-block tail
        self.canvas.create_text(450, 195, text="Custom G-block Tail:", font=("Helvetica", 12), fill="white")
        gblock_entry = ttk.Entry(self, width=80)
        gblock_entry.insert(0, self.custom_gblock_tail if self.custom_gblock_tail else "")
        self.canvas.create_window(450, 220, window=gblock_entry)

        # BSJ distance bounds
        self.canvas.create_text(450, 260, text="CircRNA Distance Filtering (Min / Max):", font=("Helvetica", 12), fill="white")
        dist_frame = tk.Frame(self, bg="#2E2E2E")
        min_entry = ttk.Entry(dist_frame, textvariable=self.min_dist, width=10)
        min_entry.pack(side="left", padx=5)
        max_entry = ttk.Entry(dist_frame, textvariable=self.max_dist, width=10)
        max_entry.pack(side="left", padx=5)
        self.canvas.create_window(450, 285, window=dist_frame)

        # Off‑target top N & mismatch cutoff
        self.canvas.create_text(450, 320, text="Off‑target Top N & Mismatch Cutoff:", font=("Helvetica", 12), fill="white")
        cutoff_frame = tk.Frame(self, bg="#2E2E2E")
        ttk.Label(cutoff_frame, text="Top N:").pack(side="left")
        ttk.Entry(cutoff_frame, textvariable=self.top_n, width=6).pack(side="left", padx=5)
        ttk.Label(cutoff_frame, text="MM cutoff:").pack(side="left", padx=(15,0))
        ttk.Entry(cutoff_frame, textvariable=self.mm_cutoff, width=6).pack(side="left", padx=5)
        self.canvas.create_window(450, 345, window=cutoff_frame)

        # Auto export
        export_check = tk.Checkbutton(self, text="Auto-Export Results to CSV", variable=self.auto_export,
                                      font=("Helvetica", 12), bg="white", fg="black")
        self.canvas.create_window(450, 390, window=export_check)

        # Buttons
        def on_confirm():
            raw_terms = motif_entry.get().strip().upper()
            self.custom_terms = [t.strip() for t in raw_terms.split(",") if t.strip()]
            self.custom_gblock_tail = gblock_entry.get().strip().upper()
            self.load_fd_scores()  # reload FD if path changed
            messagebox.showinfo("Settings Saved", "Custom settings applied.")

        def on_reset():
            self.custom_terms = []
            self.custom_gblock_tail = ""
            motif_entry.delete(0, tk.END)
            gblock_entry.delete(0, tk.END)
            self.min_dist.set(-7)
            self.max_dist.set(17)
            self.top_n.set(1000)
            self.mm_cutoff.set(5)
            self.auto_export.set(False)
            self.fd_csv_path.set("00_inputs/Raw_PFS_FD.csv")
            self.load_fd_scores()
            messagebox.showinfo("Settings Reset", "Reverted to default settings.")

        self.canvas.create_window(370, 440, window=ttk.Button(self, text="Confirm", command=on_confirm))
        self.canvas.create_window(530, 440, window=ttk.Button(self, text="Reset to Default", command=on_reset))
        self.canvas.create_window(450, 490, window=ttk.Button(self, text="Back to Menu", command=self.main_menu))

    def browse_fd_csv(self, entry_widget):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            entry_widget.delete(0, tk.END)
            entry_widget.insert(0, path)

    # ---------- Designer modes (circular / linear) ----------
    # (These are largely unchanged except they now call the backend with fd_map and no pfs_method)

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
        self.canvas.create_text(450, 50, text="Enter Upstream of Upstream (Linear Overhang)", font=("Helvetica", 18), fill="white")
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
        self.canvas.create_text(450, 50, text="Would you like to add internal exons?", font=("Helvetica", 18), fill="white")
        self.canvas.create_window(370, 100, window=ttk.Button(self, text="Yes",
                                                              command=lambda: self.get_internal_exon_loop(upstream, downstream, linearoverhang, terms)))
        self.canvas.create_window(530, 100, window=ttk.Button(self, text="No",
                                                              command=lambda: self.run_circ_analysis(upstream, downstream, linearoverhang, terms)))

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
                                                              command=lambda: self.run_circ_analysis(upstream, downstream, linearoverhang, terms)))
        self.canvas.create_window(530, 100, window=ttk.Button(self, text="No",
                                                              command=lambda: self.get_internal_exon_loop(upstream, downstream, linearoverhang, terms)))

    def run_circ_analysis(self, upstream, downstream, linearoverhang, terms):
        term_list = terms if terms else (self.custom_terms if self.custom_terms else DEFAULT_TERMS)
        gblock_tail = self.custom_gblock_tail if self.custom_gblock_tail else DEFAULT_GBLOCK_TAIL

        result = analyse_sequence(
            upstream, downstream, term_list, linearoverhang,
            gblock_tail=gblock_tail,
            fd_map=self.fd_map,
            min_dist=self.min_dist.get(),
            max_dist=self.max_dist.get()
        )
        df = result[-1]
        self.show_results_popup("Circular RNA Results", df, is_circular=True)

        if self.auto_export.get() and not df.empty:
            self.export_dataframe(df)

        if self.internal_exon_data:
            linear_df = analyse_linear_sequence(
                self.internal_exon_data, term_list, gblock_tail,
                fd_map=self.fd_map
            )
            self.show_results_popup("Internal Exon Linear Analysis", linear_df, is_circular=False)
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

            df = analyse_linear_sequence(sequence, term_list, gblock_tail, fd_map=self.fd_map)
            self.show_results_popup("Linear Analysis Results", df, is_circular=False)

            if self.auto_export.get() and not df.empty:
                self.export_dataframe(df)

            self.return_to_menu()

        self.canvas.create_window(400, 150, window=ttk.Button(self, text="Run Analysis", command=on_submit_linear))
        self.canvas.create_window(500, 150, window=ttk.Button(self, text="Cancel", command=self.main_menu))

    # ---------- Batch modes ----------
    def linear_batch_mode(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV or Excel files", "*.csv *.xlsx")])
        if not file_path:
            return
        try:
            df = pd.read_csv(file_path) if file_path.endswith('.csv') else pd.read_excel(file_path)
            if "sequence" not in df.columns:
                messagebox.showerror("Error", "Missing required 'sequence' column.")
                return
            term_list = self.custom_terms if self.custom_terms else DEFAULT_TERMS
            gblock_tail = self.custom_gblock_tail if self.custom_gblock_tail else DEFAULT_GBLOCK_TAIL
            all_dfs = []
            for i, row in df.iterrows():
                name = row.get("name", f"Entry {i+1}")
                seq = row["sequence"]
                res_df = analyse_linear_sequence(seq, term_list, gblock_tail, fd_map=self.fd_map)
                if not res_df.empty:
                    res_df.insert(0, "Sequence_Name", name)
                    all_dfs.append(res_df)
            if all_dfs:
                combined = pd.concat(all_dfs, ignore_index=True)
                self.show_results_popup("Linear Batch Results", combined, is_circular=False)
            else:
                messagebox.showinfo("No hits", "No PFS motifs found in any sequence.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to process file:\n{e}")

    def batch_circular_mode(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV or Excel files", "*.csv *.xlsx")])
        if not file_path:
            return
        try:
            df = pd.read_csv(file_path) if file_path.endswith('.csv') else pd.read_excel(file_path)
            required = ['name', 'upstream', 'downstream', 'linear_overhang']
            if not all(c in df.columns for c in required):
                messagebox.showerror("Error", f"File must contain columns: {', '.join(required)}")
                return
            term_list = self.custom_terms if self.custom_terms else DEFAULT_TERMS
            gblock_tail = self.custom_gblock_tail if self.custom_gblock_tail else DEFAULT_GBLOCK_TAIL
            all_dfs = []
            for _, row in df.iterrows():
                name = str(row['name'])
                upstream = str(row['upstream']).upper()
                downstream = str(row['downstream']).upper()
                linear_overhang = str(row['linear_overhang']).upper()
                internal_exons = str(row.get('internal_exons', '')).upper()
                _, _, _, _, circ_df = analyse_sequence(
                    upstream, downstream, term_list, linear_overhang,
                    gblock_tail=gblock_tail, fd_map=self.fd_map,
                    min_dist=self.min_dist.get(), max_dist=self.max_dist.get()
                )
                if not circ_df.empty:
                    circ_df.insert(0, "Type", "Circular")
                    circ_df.insert(0, "Sequence_Name", name)
                    all_dfs.append(circ_df)
                if internal_exons and internal_exons != 'nan':
                    lin_df = analyse_linear_sequence(internal_exons, term_list, gblock_tail, fd_map=self.fd_map)
                    if not lin_df.empty:
                        lin_df.insert(0, "Type", "Internal Exon (Linear)")
                        lin_df.insert(0, "Sequence_Name", name)
                        all_dfs.append(lin_df)
            if all_dfs:
                combined = pd.concat(all_dfs, ignore_index=True)
                self.show_results_popup("Batch Circular Results", combined, is_circular=True)
            else:
                messagebox.showinfo("No hits", "No PFS motifs found.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to process file:\n{e}")

    # ---------- Result display with off‑target button ----------
    def show_results_popup(self, title, df, is_circular=False):
        popup = tk.Toplevel(self)
        popup.title(title)
        popup.geometry("1000x500")
        frame = ttk.Frame(popup)
        frame.pack(fill="both", expand=True, padx=5, pady=5)

        if df.empty:
            ttk.Label(frame, text="No matching motifs found.").pack()
            return

        # Build Treeview
        cols = list(df.columns)
        tree = ttk.Treeview(frame, columns=cols, show="headings", height=15)
        for col in cols:
            tree.heading(col, text=col)
            tree.column(col, width=100, anchor="center")
        for _, row in df.iterrows():
            tree.insert("", "end", values=list(row))

        # Scrollbars
        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        # Button frame
        btn_frame = ttk.Frame(popup)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="Export Results", command=lambda: self.export_dataframe(df)).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Check Off‑Targets for Selected", command=lambda: self.scan_selected(tree, df, is_circular)).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Close", command=popup.destroy).pack(side="left", padx=5)

    def scan_selected(self, tree, df, is_circular):
        selection = tree.selection()
        if not selection:
            messagebox.showwarning("No selection", "Please select a row in the table.")
            return
        idx = tree.index(selection[0])  # first selected row index
        row = df.iloc[idx]
        rc_target = row["RC Target Region"]
        # Determine if it needs truncation (linear mode gives 28 nt)
        if len(rc_target) > 24:
            guide = get_24nt_guide(rc_target)
        else:
            guide = rc_target
        self.run_off_target_scan(guide)

    def run_off_target_scan(self, guide):
        # Run in thread to keep GUI responsive
        def task():
            self.progress_popup = tk.Toplevel(self)
            self.progress_popup.title("Off‑Target Scan")
            ttk.Label(self.progress_popup, text="Scanning off‑targets, please wait...").pack(padx=20, pady=10)
            progress_var = tk.StringVar(value="Starting...")
            ttk.Label(self.progress_popup, textvariable=progress_var).pack(padx=20, pady=5)
            self.progress_popup.update()

            def update_status(msg):
                progress_var.set(msg)
                self.progress_popup.update()

            try:
                results = run_off_target_for_guide(
                    guide,
                    workdir=os.getcwd(),
                    top_n=self.top_n.get(),
                    mm_cutoff=self.mm_cutoff.get(),
                    progress_callback=update_status
                )
                self.progress_popup.destroy()
                # Display results
                self.show_off_target_results(guide, results)
            except Exception as e:
                self.progress_popup.destroy()
                messagebox.showerror("Error", f"Off‑target scan failed:\n{e}")

        threading.Thread(target=task, daemon=True).start()

    def show_off_target_results(self, guide, results):
        popup = tk.Toplevel(self)
        popup.title(f"Off‑Target Hits for {guide}")
        popup.geometry("1100x500")
        if not results:
            ttk.Label(popup, text="No off‑target hits found.").pack()
            return
        # Convert list of dicts to DataFrame for easy display
        df = pd.DataFrame(results)
        # Rearrange columns for readability
        preferred_cols = [
            "gene_name", "rna_type", "pfs_5mer", "total_mismatches",
            "mismatch_intolerant", "mismatch_tolerant", "FD_score",
            "protospacer_total_hits", "protospacer_unique_records",
            "protospacer_unique_genes", "mismatch_positions_1based_all"
        ]
        available = [c for c in preferred_cols if c in df.columns]
        df = df[available + [c for c in df.columns if c not in available]]

        tree = ttk.Treeview(popup, columns=list(df.columns), show="headings", height=15)
        for col in df.columns:
            tree.heading(col, text=col)
            tree.column(col, width=100, anchor="center")
        for _, row in df.iterrows():
            tree.insert("", "end", values=list(row))
        tree.pack(fill="both", expand=True, padx=5, pady=5)

        btn_frame = ttk.Frame(popup)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="Export Results", command=lambda: self.export_dataframe(df)).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Close", command=popup.destroy).pack(side="left", padx=5)

    # ---------- Off‑target window (manual entry) ----------
    def off_target_window(self):
        win = tk.Toplevel(self)
        win.title("Off‑Target Scan (Manual)")
        win.geometry("400x200")
        ttk.Label(win, text="Enter 24‑nt guide sequence:").pack(pady=10)
        entry = ttk.Entry(win, width=40)
        entry.pack(pady=5)
        def on_scan():
            guide = entry.get().strip().upper().replace("U", "T")
            if len(guide) != 24:
                messagebox.showwarning("Invalid", "Guide must be 24 nucleotides.")
                return
            win.destroy()
            self.run_off_target_scan(guide)
        ttk.Button(win, text="Scan", command=on_scan).pack(pady=10)

    # ---------- Export helper ----------
    def export_dataframe(self, df):
        if df.empty:
            messagebox.showwarning("No Data", "There is no data to export.")
            return
        file_path = filedialog.asksaveasfilename(defaultextension=".csv",
                                                 filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx")])
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

    def return_to_menu(self):
        self.internal_exon_data = ""
        self.after(1500, self.main_menu)


if __name__ == "__main__":
    app = SequenceApp()
    app.mainloop()