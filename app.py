"""
Fluorescence Analysis GUI
Wraps fluorescence_analysis.main with a Tkinter interface.
"""

import shutil
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk

from PIL import Image, ImageTk

import fluorescence_analysis as fa


# ── helpers ───────────────────────────────────────────────────────────────────

THUMB = (320, 240)


def _run_analysis(img_dir: Path, log: scrolledtext.ScrolledText,
                  btn_run: tk.Button, after_cb):
    """Run analysis in a background thread, pipe stdout to the log widget."""
    import io, sys

    class _Tee(io.TextIOBase):
        def write(self, s):
            log.after(0, lambda: (_append_log(log, s)))
            return len(s)

    old_stdout = sys.stdout
    sys.stdout = _Tee()
    try:
        fa.main(str(img_dir))
    except SystemExit as e:
        log.after(0, lambda: _append_log(log, f"\nERROR: {e}\n"))
    except Exception as e:
        log.after(0, lambda: _append_log(log, f"\nERROR: {e}\n"))
    finally:
        sys.stdout = old_stdout
        log.after(0, lambda: (btn_run.config(state="normal"), after_cb()))


def _append_log(log: scrolledtext.ScrolledText, text: str):
    log.config(state="normal")
    log.insert("end", text)
    log.see("end")
    log.config(state="disabled")


# ── gallery ───────────────────────────────────────────────────────────────────

def _build_gallery(frame: tk.Frame, png_paths: list[Path]):
    for w in frame.winfo_children():
        w.destroy()

    cols = 4
    for idx, p in enumerate(png_paths):
        img = Image.open(p)
        img.thumbnail(THUMB)
        photo = ImageTk.PhotoImage(img)

        cell = tk.Frame(frame, bd=1, relief="solid")
        cell.grid(row=idx // cols, column=idx % cols, padx=4, pady=4)

        lbl = tk.Label(cell, image=photo, cursor="hand2")
        lbl.image = photo          # keep reference
        lbl.pack()
        tk.Label(cell, text=p.name, font=("Arial", 7), wraplength=THUMB[0]).pack()

        lbl.bind("<Button-1>", lambda e, path=p: _enlarge(path))


def _enlarge(path: Path):
    top = tk.Toplevel()
    top.title(path.name)
    img = Image.open(path)
    img.thumbnail((900, 700))
    photo = ImageTk.PhotoImage(img)
    tk.Label(top, image=photo).pack()
    top.photo = photo              # keep reference


# ── instructions panel ───────────────────────────────────────────────────────

_INSTRUCTIONS = [
    ("Overview", None, [
        "Channel-agnostic fluorescence quantifier — works with any combination of",
        "GFP, RFP, DAPI, CY5, mCherry, YFP, CFP, TRITC, FITC, Alexa dyes, and more.",
        "No fixed pairing required: bring 1, 2, 3+ channels and any number of replicates.",
        "Computes ImageJ-equivalent metrics (Mean, IntDen) and all pairwise channel ratios.",
    ]),
    ("File Format", None, [
        "• Extension : .tif  (TIFF only — .tiff / .png / .jpg are NOT supported)",
        "• Bit depth : 8-bit or 16-bit grayscale, or RGB / RGBA TIFF",
        "• No fixed channel pairing needed — any subset of channels per well is fine",
    ]),
    ("Recognised Channel Tags", None, [
        "The channel tag must appear in the filename, preceded by an underscore.",
        "Matching is case-insensitive.",
        "",
        "  GFP   RFP   DAPI   CY5   MCHERRY   YFP   CFP",
        "  TRITC   FITC   A488   A555   A647   TXR",
        "",
        "Any file whose stem does NOT contain one of these tags is skipped with a warning.",
        "To add a custom tag, edit _CHANNEL_RE in fluorescence_analysis.py.",
    ]),
    ("File Naming Convention", None, [
        "Pattern:  <WellName>_<CHANNEL>[_<replicate>].tif",
        "",
        "Rules:",
        "  • WellName can contain underscores (treatment, inhibitor, etc.)",
        "  • Channel tag must directly follow an underscore",
        "  • Optional replicate number goes at the very end:  _1  _2  _3 …",
        "  • Multiple replicates of the same channel are averaged before quantification",
        "  • A well can have only one channel, or many — all are processed",
        "",
        "Valid examples:",
        "  MAPK10_GFP_1.tif            single channel, replicate 1",
        "  MAPK10_GFP_2.tif            same well, replicate 2  (averaged with rep 1)",
        "  MAPK10_RFP_1.tif            second channel for the same well",
        "  MAPK10_DAPI_1.tif           third channel — DAPI added freely",
        "  MAPK10_JNK_GFP_1.tif       inhibitor JNK, GFP channel, replicate 1",
        "  PCAG_CY5.tif               no replicate suffix — treated as replicate 1",
        "  SAMPLE_A488_1.tif          Alexa-488 channel",
    ]),
    ("Well / Condition Naming", None, [
        "Pattern:  <TreatmentGroup>[_<Inhibitor>]_<Channel>[_<Replicate>]",
        "",
        "Recognised inhibitor keywords (case-insensitive, must follow an underscore):",
        "  WORTMAN  |  UO126  |  SB90  |  JNK",
        "",
        "Examples:",
        "  MAPK10_GFP_1.tif              → treatment: MAPK10,          inhibitor: Control",
        "  MAPK10_JNK_GFP_2.tif         → treatment: MAPK10,          inhibitor: JNK",
        "  MAPK10_GENTECIN_GFP_1.tif    → treatment: MAPK10_GENTECIN, inhibitor: Control",
        "  MAPK10_GENTECIN_UO126_RFP_1  → treatment: MAPK10_GENTECIN, inhibitor: UO126",
    ]),
    ("Folder Structure", None, [
        "Place ALL .tif files flat inside one folder (no sub-folders):",
        "",
        "  my_experiment/",
        "  ├── MAPK10_GFP_1.tif",
        "  ├── MAPK10_GFP_2.tif          ← replicate 2, averaged with rep 1",
        "  ├── MAPK10_RFP_1.tif",
        "  ├── MAPK10_DAPI_1.tif         ← extra channel, no pairing needed",
        "  ├── MAPK10_JNK_GFP_1.tif",
        "  ├── PCAG_CY5.tif              ← single-channel well, no problem",
        "  └── ...",
        "",
        "Results are written to:  <parent_of_selected_folder>/results/",
    ]),
    ("Outputs", None, [
        "CSV files:",
        "  fluorescence_results.csv   — per-well: Mean + IntDen for every channel,",
        "                               plus all pairwise ratios (e.g. GFP_RFP_ratio)",
        "  condition_summary.csv      — mean ± SD per condition group",
        "",
        "PNG plots (per well):",
        "  <WellName>_analysis.png    — channel images, histograms, pairwise scatter",
        "",
        "PNG plots (summary):",
        "  summary_intden.png              — all channels side-by-side per well",
        "  summary_<A>_<B>_ratio.png       — one plot per channel pair",
        "  condition_<channel>_intden.png  — per-channel condition comparison",
        "  condition_<A>_<B>_ratio.png     — per-ratio condition comparison",
        "  condition_heatmap.png           — z-scored heatmap across all metrics",
        "  replicate_consistency.png       — rep1 vs rep2 scatter",
        "  inhibitor_comparison_<group>.png",
    ]),
    ("Quick Start", None, [
        "1. Click  Browse…  and select the folder containing your .tif files.",
        "2. Click  ▶ Run Analysis  — progress appears in the Progress Log tab.",
        "3. When complete, the PNG Gallery tab opens automatically.",
        "4. Click  ⬇ Download CSV  to save the results spreadsheet anywhere.",
    ]),
]


def _build_instructions(parent: tk.Frame):
    canvas = tk.Canvas(parent, highlightthickness=0)
    vsb = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    inner = tk.Frame(canvas, padx=18, pady=12)
    win = canvas.create_window((0, 0), window=inner, anchor="nw")

    def _on_frame(e):  canvas.configure(scrollregion=canvas.bbox("all"))
    def _on_canvas(e): canvas.itemconfig(win, width=e.width)
    inner.bind("<Configure>", _on_frame)
    canvas.bind("<Configure>", _on_canvas)

    # bind mousewheel
    def _scroll(e): canvas.yview_scroll(int(-1 * (e.delta / 120)), "units")
    canvas.bind_all("<MouseWheel>", _scroll)

    tk.Label(inner, text="Fluorescence Analysis — User Guide",
             font=("Arial", 14, "bold"), fg="#2ca02c").pack(anchor="w", pady=(0, 8))

    for title, _bg, lines in _INSTRUCTIONS:
        # section header
        tk.Label(inner, text=title, font=("Arial", 10, "bold"),
                 fg="#2ca02c").pack(anchor="w", pady=(10, 2))
        tk.Frame(inner, height=1, bg="#cccccc").pack(fill="x", pady=(0, 4))

        block_text = "\n".join(lines)
        tk.Label(inner, text=block_text, font=("Courier", 9),
                 justify="left", anchor="w", wraplength=700,
                 bg="#f5f5f5", relief="flat", padx=8, pady=6).pack(
                     fill="x", anchor="w")


# ── main window ───────────────────────────────────────────────────────────────

def build_ui():
    root = tk.Tk()
    root.title("Fluorescence Analysis")
    root.resizable(True, True)
    root.state("zoomed")   # fullscreen windowed (maximised, keeps title bar)

    state = {"img_dir": None, "out_dir": None}

    # ── top bar ───────────────────────────────────────────────────────────────
    top = tk.Frame(root, pady=6, padx=8)
    top.pack(fill="x")

    dir_var = tk.StringVar(value="No folder selected")
    tk.Label(top, text="Image folder:").pack(side="left")
    tk.Label(top, textvariable=dir_var, relief="sunken", width=50,
             anchor="w").pack(side="left", padx=4)

    def pick_folder():
        d = filedialog.askdirectory(title="Select image folder")
        if d:
            state["img_dir"] = Path(d)
            dir_var.set(d)
            btn_run.config(state="normal")

    tk.Button(top, text="Browse…", command=pick_folder).pack(side="left")

    # ── notebook ──────────────────────────────────────────────────────────────
    nb = ttk.Notebook(root)
    nb.pack(fill="both", expand=True, padx=8, pady=4)

    # Tab 0 – instructions
    instr_frame = tk.Frame(nb)
    nb.add(instr_frame, text="📋 Instructions")
    _build_instructions(instr_frame)

    # Tab 1 – log
    log_frame = tk.Frame(nb)
    nb.add(log_frame, text="Progress Log")
    log = scrolledtext.ScrolledText(log_frame, state="disabled",
                                    font=("Courier", 9), wrap="word")
    log.pack(fill="both", expand=True)

    # Tab 3 – gallery
    gallery_outer = tk.Frame(nb)
    nb.add(gallery_outer, text="PNG Gallery")

    canvas = tk.Canvas(gallery_outer)
    vsb = ttk.Scrollbar(gallery_outer, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    gallery_frame = tk.Frame(canvas)
    canvas_win = canvas.create_window((0, 0), window=gallery_frame, anchor="nw")

    def _on_frame_configure(e):
        canvas.configure(scrollregion=canvas.bbox("all"))

    def _on_canvas_configure(e):
        canvas.itemconfig(canvas_win, width=e.width)

    gallery_frame.bind("<Configure>", _on_frame_configure)
    canvas.bind("<Configure>", _on_canvas_configure)

    # ── bottom bar ────────────────────────────────────────────────────────────
    bot = tk.Frame(root, pady=6, padx=8)
    bot.pack(fill="x")

    def download_csv():
        src = state.get("out_dir")
        if not src:
            messagebox.showwarning("No results", "Run analysis first.")
            return
        csv_src = src / "fluorescence_results.csv"
        if not csv_src.exists():
            messagebox.showwarning("Not found", f"{csv_src} not found.")
            return
        dest = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile="fluorescence_results.csv",
        )
        if dest:
            shutil.copy(csv_src, dest)
            messagebox.showinfo("Saved", f"CSV saved to:\n{dest}")

    btn_csv = tk.Button(bot, text="⬇ Download CSV", command=download_csv,
                        state="disabled")
    btn_csv.pack(side="right", padx=4)

    def after_analysis():
        out_dir = state["out_dir"]
        pngs = sorted(out_dir.glob("*.png"))
        _build_gallery(gallery_frame, pngs)
        nb.select(2)
        btn_csv.config(state="normal")

    def run_analysis():
        if not state["img_dir"]:
            return
        state["out_dir"] = state["img_dir"].parent / "results"
        btn_run.config(state="disabled")
        log.config(state="normal"); log.delete("1.0", "end"); log.config(state="disabled")
        t = threading.Thread(
            target=_run_analysis,
            args=(state["img_dir"], log, btn_run, after_analysis),
            daemon=True,
        )
        t.start()
        nb.select(1)

    btn_run = tk.Button(bot, text="▶ Run Analysis", command=run_analysis,
                        state="disabled", bg="#2ca02c", fg="white",
                        font=("Arial", 10, "bold"))
    btn_run.pack(side="left", padx=4)

    root.mainloop()


if __name__ == "__main__":
    build_ui()
