from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import streamlit as st


REPO_ROOT = Path(__file__).resolve().parent
R_SCRIPT = REPO_ROOT / "glioma_analysis.R"
DEFAULT_OUTPUT = REPO_ROOT / "analysis_outputs"


st.set_page_config(
    page_title="Glioma Multi-Omics Analysis Suite",
    page_icon="🧬",
    layout="wide",
)

st.markdown(
    """
    <style>
      :root {
        --ink: #1c2433;
        --muted: #5a6a7a;
        --accent: #1f4e5f;
        --accent-strong: #0f766e;
        --paper: #f7f3ea;
        --panel: #ffffff;
        --stroke: #e1e6ee;
      }

      .stApp {
        background: radial-gradient(1200px 800px at 15% 10%, #f7efe0 0%, transparent 60%),
                    radial-gradient(900px 700px at 85% 15%, #eef3f6 0%, transparent 55%),
                    linear-gradient(180deg, #f7f5f0 0%, #f1f4f8 60%, #f5f2ed 100%);
        font-family: "Iowan Old Style", "Palatino Linotype", "Book Antiqua", "Garamond", serif;
        color: var(--ink);
      }

      h1, h2, h3, h4 {
        font-family: "Iowan Old Style", "Palatino Linotype", "Book Antiqua", "Garamond", serif;
        letter-spacing: 0.2px;
      }

      .header-block {
        background: rgba(255, 255, 255, 0.92);
        border: 1px solid var(--stroke);
        border-radius: 18px;
        padding: 22px 26px;
        box-shadow: 0 18px 40px rgba(15, 23, 42, 0.08);
        margin-bottom: 20px;
      }

      .badge {
        display: inline-block;
        padding: 6px 12px;
        border-radius: 999px;
        background: rgba(31, 78, 95, 0.12);
        color: var(--accent);
        font-weight: 600;
        letter-spacing: 0.3px;
        font-size: 12px;
        text-transform: uppercase;
      }

      .card {
        background: rgba(255, 255, 255, 0.9);
        border: 1px solid var(--stroke);
        border-radius: 16px;
        padding: 16px 18px;
        box-shadow: 0 10px 30px rgba(15, 23, 42, 0.08);
      }

      .stat {
        font-size: 28px;
        font-weight: 700;
        color: var(--accent-strong);
        margin-top: 4px;
      }

      .muted {
        color: var(--muted);
      }

      .fade-in {
        animation: fadeIn 0.8s ease-in both;
      }

      @keyframes fadeIn {
        from { opacity: 0; transform: translateY(8px); }
        to { opacity: 1; transform: translateY(0); }
      }

      /* Fix tab text visibility */
      .stTabs [data-baseweb="tab-list"] button {
        color: #1c2433 !important;
        font-weight: 600 !important;
      }
      
      .stTabs [data-baseweb="tab-list"] button[aria-selected="true"] {
        color: #0f766e !important;
        font-weight: 700 !important;
      }

      /* Make all text more visible */
      p, span, div, label {
        color: #1c2433;
      }

      /* Sidebar text */
      .stSidebar, .stSidebar label, .stSidebar p, .stSidebar span {
        color: #1c2433 !important;
      }

      /* Dark sidebar: ensure readable text */
      section[data-testid="stSidebar"] {
        color: #e5e7eb !important;
      }

      section[data-testid="stSidebar"] * {
        color: #e5e7eb !important;
      }

      section[data-testid="stSidebar"] label,
      section[data-testid="stSidebar"] p,
      section[data-testid="stSidebar"] span,
      section[data-testid="stSidebar"] div {
        color: #e5e7eb !important;
      }

      section[data-testid="stSidebar"] small,
      section[data-testid="stSidebar"] .stMarkdown,
      section[data-testid="stSidebar"] .stCaption {
        color: #cbd5f5 !important;
      }

      /* Headers */
      h1, h2, h3, h4, h5, h6 {
        color: #1c2433 !important;
      }

      /* Expander text */
      .streamlit-expanderHeader {
        color: #1c2433 !important;
        font-weight: 600 !important;
      }

      /* Download buttons: improve contrast */
      .stDownloadButton button {
        background: linear-gradient(180deg, #1f2937 0%, #111827 100%) !important;
        color: #f8fafc !important;
        border: 1px solid #334155 !important;
        box-shadow: 0 10px 20px rgba(15, 23, 42, 0.18) !important;
        font-weight: 600 !important;
        letter-spacing: 0.15px !important;
      }

      .stDownloadButton button:hover {
        background: linear-gradient(180deg, #1f4e5f 0%, #0f3b4a 100%) !important;
        color: #ffffff !important;
        border-color: #1f4e5f !important;
        box-shadow: 0 12px 24px rgba(15, 23, 42, 0.22) !important;
      }

      .stDownloadButton button:active {
        transform: translateY(0);
        box-shadow: 0 6px 12px rgba(15, 23, 42, 0.18) !important;
      }

      .stDownloadButton button:disabled {
        background: #2b3441 !important;
        color: #9ca3af !important;
        border-color: #3b4250 !important;
        box-shadow: none !important;
      }

      .stDownloadButton button span,
      .stDownloadButton button p,
      .stDownloadButton button div {
        color: #f8fafc !important;
      }

      .stDownloadButton button:disabled span,
      .stDownloadButton button:disabled p,
      .stDownloadButton button:disabled div {
        color: #9ca3af !important;
      }

      /* Tooltip/help popovers */
      [role="tooltip"],
      .stTooltip {
        color: #f8fafc !important;
        background: #111827 !important;
        border: 1px solid #1f2937 !important;
        box-shadow: 0 10px 24px rgba(15, 23, 42, 0.2) !important;
      }

      [role="tooltip"] *,
      .stTooltip * {
        color: #f8fafc !important;
      }

      /* BaseWeb tooltip container used by Streamlit */
      [data-baseweb="tooltip"] {
        color: #f8fafc !important;
        background: #111827 !important;
        border: 1px solid #1f2937 !important;
        box-shadow: 0 10px 24px rgba(15, 23, 42, 0.2) !important;
      }

      [data-baseweb="tooltip"] * {
        color: #f8fafc !important;
      }
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <div class="header-block fade-in">
      <span class="badge">Glioma Research Toolkit</span>
      <h1>Glioma Molecular Profiling and Prediction Suite</h1>
      <p class="muted">
        A streamlined academic interface for DESeq2, pathway enrichment, GSVA,
        and multi-model classification benchmarking with hierarchical (TYPE + GRADE) prediction.
      </p>
    </div>
    """,
    unsafe_allow_html=True,
)

# Model information expander
with st.expander("ℹ️ About the Enhanced Classification Pipeline"):
    st.markdown(
        """
        ### Two-Stage Hierarchical Classification
        
        **Stage 1 - TYPE Classification:**
        - Predicts glioma type: Astrocytoma, Oligodendroglioma, Oligodendroastrocytoma, Glioblastoma
        
        **Stage 2 - GRADE Classification:**
        - Predicts tumor grade: Grade II, Grade III, Grade IV (Primary/Secondary), Recurrent
        
        ### Enhanced ML Pipeline Features
        
        | Feature | Description |
        |---------|-------------|
        | **Feature Selection** | Combines variance, correlation, and ANOVA F-statistics |
        | **Class Balancing** | Oversampling of minority classes with noise injection |
        | **Hyperparameter Tuning** | Grid search for optimal model parameters |
        | **8 Models Tested** | Random Forest, XGBoost, SVM, KNN, Neural Network, Elastic Net, LDA, Ensemble |
        | **Ensemble Voting** | Majority voting across all models for robust predictions |
        """
    )

with st.sidebar:
    st.markdown("### Run Configuration")
    counts_file = st.file_uploader("Counts CSV", type=["csv"])
    meta_file = st.file_uploader("Metadata CSV", type=["csv"])
    output_dir = st.text_input("Output directory", value=str(DEFAULT_OUTPUT))
    run_visuals = st.checkbox("Generate extended visualizations", value=False)
    run_button = st.button("Run pipeline", type="primary", use_container_width=True)

st.markdown("### Project Context")
st.markdown(
    """
    <div class="card">
      <p class="muted">
        Upload the raw counts matrix and metadata table, then run the pipeline.
        Outputs are saved locally in the specified output directory, including
        DESeq2 results, model benchmarks, and summary plots.
      </p>
    </div>
    """,
    unsafe_allow_html=True,
)


def _write_upload(upload, destination: Path) -> Path:
    destination.write_bytes(upload.getbuffer())
    return destination


def _check_rscript() -> bool:
    """Check if Rscript is available in the system PATH."""
    import shutil
    return shutil.which("Rscript") is not None


def _run_pipeline(counts_path: Path, meta_path: Path, out_dir: Path, make_viz: bool) -> tuple[int, str]:
    if not _check_rscript():
        return 1, "Error: Rscript not found. Please install R from https://cran.r-project.org/ or via 'brew install r'"
    
    cmd = [
        "Rscript",
        str(R_SCRIPT),
        "--counts",
        str(counts_path),
        "--meta",
        str(meta_path),
        "--out",
        str(out_dir),
    ]
    if make_viz:
        cmd.append("--viz")
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, (result.stdout + "\n" + result.stderr).strip()


if run_button:
    if counts_file is None or meta_file is None:
        st.error("Please upload both the counts CSV and metadata CSV.")
    else:
        output_path = Path(output_dir).expanduser()
        input_dir = output_path / "inputs"
        input_dir.mkdir(parents=True, exist_ok=True)

        counts_path = _write_upload(counts_file, input_dir / f"counts_{counts_file.name}")
        meta_path = _write_upload(meta_file, input_dir / f"meta_{meta_file.name}")

        with st.spinner("Running analysis pipeline. This may take a few minutes."):
            code, log_output = _run_pipeline(counts_path, meta_path, output_path, run_visuals)

        if code != 0:
            st.error("Pipeline failed. Review the log output for details.")
        else:
            st.success("Pipeline completed successfully.")

        with st.expander("Pipeline log"):
            st.code(log_output or "No logs captured.")

st.markdown("### Results Overview")

results_dir = Path(output_dir).expanduser() / "analysis_results"
accuracy_path = results_dir / "model_accuracy.csv"
type_accuracy_path = results_dir / "model_accuracy_TYPE.csv"
grade_accuracy_path = results_dir / "model_accuracy_GRADE.csv"
hierarchical_path = results_dir / "hierarchical_comparison.csv"
best_model_path = results_dir / "best_model.txt"

if accuracy_path.exists():
    accuracy_df = pd.read_csv(accuracy_path)
    best_model = None
    if best_model_path.exists():
        best_model = best_model_path.read_text().strip().replace("BestModel=", "")

    # Load TYPE and GRADE specific results
    type_df = pd.read_csv(type_accuracy_path) if type_accuracy_path.exists() else None
    grade_df = pd.read_csv(grade_accuracy_path) if grade_accuracy_path.exists() else None
    hierarchical_df = pd.read_csv(hierarchical_path) if hierarchical_path.exists() else None

    # Get best accuracies for TYPE and GRADE
    type_best_acc = type_df["Accuracy"].max() if type_df is not None and "Accuracy" in type_df else None
    grade_best_acc = grade_df["Accuracy"].max() if grade_df is not None and "Accuracy" in grade_df else None
    overall_best_acc = accuracy_df["Accuracy"].max() if "Accuracy" in accuracy_df else None
    
    # Get best model names for each stage
    type_best_model = type_df.loc[type_df["Accuracy"].idxmax(), "Model"] if type_df is not None and "Accuracy" in type_df else "N/A"
    grade_best_model = grade_df.loc[grade_df["Accuracy"].idxmax(), "Model"] if grade_df is not None and "Accuracy" in grade_df else "N/A"

    # Header metrics - Summary Cards
    st.markdown("#### 🎯 Classification Performance Summary")
    
    metric_cols = st.columns(4)
    with metric_cols[0]:
        st.markdown(
            f"""
            <div class="card fade-in">
              <div class="muted">🏆 Best Overall Model</div>
              <div class="stat" style="font-size: 16px;">{best_model or "Available"}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    with metric_cols[1]:
        type_acc_display = f"{type_best_acc:.1%}" if type_best_acc is not None else "n/a"
        st.markdown(
            f"""
            <div class="card fade-in">
              <div class="muted">🧬 TYPE Accuracy</div>
              <div class="stat">{type_acc_display}</div>
              <div class="muted" style="font-size: 11px;">{type_best_model}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    with metric_cols[2]:
        grade_acc_display = f"{grade_best_acc:.1%}" if grade_best_acc is not None else "n/a"
        st.markdown(
            f"""
            <div class="card fade-in">
              <div class="muted">📊 GRADE Accuracy</div>
              <div class="stat">{grade_acc_display}</div>
              <div class="muted" style="font-size: 11px;">{grade_best_model}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    with metric_cols[3]:
        overall_acc_display = f"{overall_best_acc:.1%}" if overall_best_acc is not None else "n/a"
        num_models = len(accuracy_df) if accuracy_df is not None else 0
        st.markdown(
            f"""
            <div class="card fade-in">
              <div class="muted">📈 Best Accuracy</div>
              <div class="stat">{overall_acc_display}</div>
              <div class="muted" style="font-size: 11px;">{num_models} models tested</div>
            </div>
            """,
            unsafe_allow_html=True,
        )

    # Hierarchical Comparison - Show first if available
    if hierarchical_df is not None:
        st.markdown("#### 🔄 Hierarchical vs Direct Classification Comparison")
        st.markdown(
            """
            <div class="card">
              <p class="muted">
                Comparing the two-stage hierarchical approach (TYPE → GRADE) against direct full classification.
              </p>
            </div>
            """,
            unsafe_allow_html=True,
        )
        # Highlight best approach
        st.dataframe(
            hierarchical_df.style.highlight_max(subset=["Accuracy"], color="#d4edda"),
            use_container_width=True, 
            hide_index=True
        )

    # Create tabs for TYPE and GRADE results
    tab1, tab2, tab3 = st.tabs(["🧬 Stage 1: TYPE", "📊 Stage 2: GRADE", "📋 All Models"])
    
    with tab1:
        st.markdown("#### Glioma TYPE Classification")
        st.markdown(
            """
            <div class="card">
              <p class="muted">
                Classifies tumors into: <strong>Astrocytoma</strong>, <strong>Oligodendroglioma</strong>, 
                <strong>Oligodendroastrocytoma</strong>, or <strong>Glioblastoma</strong>
              </p>
            </div>
            """,
            unsafe_allow_html=True,
        )
        if type_df is not None:
            # Sort by accuracy descending
            type_df_sorted = type_df.sort_values("Accuracy", ascending=False)
            st.dataframe(
                type_df_sorted.style.highlight_max(subset=["Accuracy"], color="#d4edda")
                                    .format({"Accuracy": "{:.4f}"}),
                use_container_width=True, 
                hide_index=True
            )
            
            # Show accuracy chart
            st.bar_chart(type_df_sorted.set_index("Model")["Accuracy"])

    with tab2:
        st.markdown("#### Glioma GRADE Classification")
        st.markdown(
            """
            <div class="card">
              <p class="muted">
                Classifies tumors into: <strong>Grade II</strong>, <strong>Grade III</strong>, 
                <strong>Grade IV (Primary/Secondary)</strong>, or <strong>Recurrent</strong>
              </p>
            </div>
            """,
            unsafe_allow_html=True,
        )
        if grade_df is not None:
            # Sort by accuracy descending
            grade_df_sorted = grade_df.sort_values("Accuracy", ascending=False)
            st.dataframe(
                grade_df_sorted.style.highlight_max(subset=["Accuracy"], color="#d4edda")
                                     .format({"Accuracy": "{:.4f}"}),
                use_container_width=True, 
                hide_index=True
            )
            
            # Show accuracy chart
            st.bar_chart(grade_df_sorted.set_index("Model")["Accuracy"])

    with tab3:
        st.markdown("#### All Models Ranked by Accuracy")
        st.markdown(
            """
            <div class="card">
              <p class="muted">
                Complete ranking of all models across TYPE, GRADE, and DIRECT classification stages.
              </p>
            </div>
            """,
            unsafe_allow_html=True,
        )
        # Sort by accuracy descending
        accuracy_df_sorted = accuracy_df.sort_values("Accuracy", ascending=False)
        st.dataframe(
            accuracy_df_sorted.style.highlight_max(subset=["Accuracy"], color="#d4edda")
                                    .format({"Accuracy": "{:.4f}"}),
            use_container_width=True, 
            hide_index=True
        )

else:
    st.info("🔬 Run the pipeline to populate results. Upload your counts and metadata CSV files in the sidebar.")

st.markdown("### 📥 Downloadable Reports")

# Create columns for download buttons
if results_dir.exists():
    dl_cols = st.columns(4)
    
    # Accuracy Summary CSV
    with dl_cols[0]:
        if accuracy_path.exists():
            st.download_button(
                "📊 All Models",
                data=accuracy_path.read_bytes(),
                file_name="model_accuracy.csv",
                mime="text/csv",
                use_container_width=True,
                help="Download accuracy results for all models"
            )
        else:
            st.caption("Not available")
    
    # TYPE Accuracy CSV
    with dl_cols[1]:
        if type_accuracy_path.exists():
            st.download_button(
                "🧬 TYPE Results",
                data=type_accuracy_path.read_bytes(),
                file_name="model_accuracy_TYPE.csv",
                mime="text/csv",
                use_container_width=True,
                help="Download TYPE classification results"
            )
        else:
            st.caption("Not available")
    
    # GRADE Accuracy CSV
    with dl_cols[2]:
        if grade_accuracy_path.exists():
            st.download_button(
                "📈 GRADE Results",
                data=grade_accuracy_path.read_bytes(),
                file_name="model_accuracy_GRADE.csv",
                mime="text/csv",
                use_container_width=True,
                help="Download GRADE classification results"
            )
        else:
            st.caption("Not available")
    
    # Analysis Summary PDF
    with dl_cols[3]:
        summary_pdf = results_dir / "analysis_summary.pdf"
        if summary_pdf.exists():
            st.download_button(
                "📄 Summary PDF",
                data=summary_pdf.read_bytes(),
                file_name=summary_pdf.name,
                mime="application/pdf",
                use_container_width=True,
                help="Download analysis summary report"
            )
        else:
            st.caption("Not available")

    # Hierarchical Comparison CSV (separate row)
    if hierarchical_path.exists():
        st.download_button(
            "📋 Download Hierarchical Comparison Report",
            data=hierarchical_path.read_bytes(),
            file_name="hierarchical_comparison.csv",
            mime="text/csv",
            help="Compare hierarchical vs direct classification approaches"
        )
        
    # Complete analysis RDS (for R users)
    complete_rds = results_dir / "complete_analysis.rds"
    if complete_rds.exists():
        st.download_button(
            "💾 Download Complete Analysis (RDS)",
            data=complete_rds.read_bytes(),
            file_name="complete_analysis.rds",
            mime="application/octet-stream",
            help="Download full analysis object for R (includes all models and data)"
        )
else:
    st.caption("📁 Reports will appear after the pipeline completes.")

# Footer
st.markdown("---")
st.markdown(
    """
    <div style="text-align: center; color: var(--muted); font-size: 12px;">
        <p>🧬 Glioma Multi-Omics Analysis Suite | Enhanced ML Pipeline with 8 Models</p>
        <p>Features: Hierarchical Classification • Hyperparameter Tuning • Ensemble Voting • Class Balancing</p>
    </div>
    """,
    unsafe_allow_html=True,
)
