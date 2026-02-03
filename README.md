# CellMaster

**CellMaster** is an intelligent automated pipeline for single-cell RNA sequencing (scRNA-seq) cell-type annotation. It leverages Large Language Models (LLMs) to iteratively refine cell-type annotations through a multi-agent system that includes hypothesis generation, marker gene proposal, visualization, and evaluation.

## 🌟 Features

- **Automated Cell-Type Annotation**: Intelligent, iterative annotation of cell types in scRNA-seq data
- **LLM-Powered Analysis**: Utilizes advanced language models (GPT-4, o1-mini) for biological reasoning
- **Multi-Agent Architecture**: Coordinated agents for hypothesis, experiment, environment, and evaluation
- **Flexible Configuration**: Support for multiple species (human, mouse) and clustering methods
- **Visualization**: Automatic generation of UMAP plots and dotplots for validation
- **Iterative Refinement**: Progressive improvement of annotations through multiple iterations
- **Marker Gene Discovery**: Intelligent proposal and validation of cell-type-specific marker genes

## 🚀 Quick Start

### Demo / Vignette

**See [`testing.ipynb`](testing.ipynb) for a complete demonstration of CellMaster's capabilities.** This Jupyter notebook provides a step-by-step example of cell-type annotation on retina and PBMC datasets.

### UI Version

For a graphical user interface version of CellMaster, check out **[@AnonymousGym/CellMaster-UI](https://github.com/AnonymousGym/CellMaster-UI)**.

## 📋 Prerequisites

- Python 3.11+
- OpenAI API key or compatible LLM API access
- Single-cell RNA-seq data in `.h5ad` format (AnnData)

## 🔧 Installation

### 1. Clone the Repository

```bash
git clone https://github.com/AnonymousGym/CellMaster.git
cd CellMaster
```

### 2. Set Up Environment

Create and activate a conda environment:

```bash
conda create -n cellmaster python=3.11
conda activate cellmaster
```

Install dependencies:

```bash
pip install -r requirements.txt
```

### 3. Configure API Key

Replace the `OPENAI_API_KEY` in `/config/settings.py` with your own API key:

```python
# config/settings.py
OPENAI_API_KEY = 'your-api-key-here'
```

### 4. Download Required Data Files

Download the example datasets from:
- [Google Drive Link](https://drive.google.com/file/d/1qT5qBbUYP3BvXNL7NY6zqEr-pp13zsdW/view?usp=sharing)

Place the downloaded files in an `uploads/` directory in the project root.

## 💻 Usage

### Basic Usage with Python API

```python
from cli import run_cell_annotation
import scanpy as sc

# Load your scRNA-seq data
adata = sc.read_h5ad("path/to/your/data.h5ad")

# Define your research hypothesis
hypothesis = """
You are a bioinformatics expert in single-cell RNA sequencing analysis.
Dataset: This is a scRNA dataset of [X] cells and [Y] genes from [tissue/cell type].
Project Description: [Describe your research goals and expected cell types]
"""

# Run cell annotation
annotated_adata = run_cell_annotation(
    input_adata=adata,
    initial_hypothesis=hypothesis,
    output_column="cell_type_annotation",
    API_key="your-api-key-here",
    model_name="gpt-4o"
)

# Access annotations
print(annotated_adata.obs["cell_type_annotation"])
```

### Advanced Configuration

```python
annotated_adata = run_cell_annotation(
    input_adata=adata,
    initial_hypothesis=hypothesis,
    output_column="cell_type_annotation",
    API_key="your-api-key-here",
    model_name="gpt-4o",
    model_provider="openai",
    original_grouping="leiden",  # clustering column name
    species="human",  # or "mouse"
    input_dir="CellMaster",
    output_dir="CellMaster",
    log_file_path="CellMaster/annotation_log.txt"
)
```

## 📂 Project Structure

```
CellMaster/
├── agents/                      # Multi-agent system components
│   ├── hypothesis_agent/        # Hypothesis generation and refinement
│   ├── experiment_agent/        # Marker gene proposal
│   ├── environment_agent/       # Visualization (dotplots, UMAP)
│   └── evaluation_agent/        # Annotation evaluation
├── config/                      # Configuration files
│   └── settings.py             # API keys and settings
├── utils/                       # Utility functions
│   ├── LLM.py                  # LLM interaction utilities
│   ├── liver_process_toolkit.py # Processing utilities
│   └── traj_util.py            # Trajectory analysis utilities
├── cli.py                       # Command-line interface
├── pipeline.py                  # Main pipeline orchestration
├── testing.ipynb               # Demo notebook (START HERE!)
├── requirements.txt            # Python dependencies
├── dependencies.yaml           # Additional dependencies
└── README.md                   # This file
```

## 🔬 How It Works

CellMaster uses a four-agent system that iteratively refines cell-type annotations:

1. **Hypothesis Agent**: Generates and refines hypotheses about cell types based on marker genes and prior knowledge
2. **Experiment Agent**: Proposes marker genes for validation of specific cell types
3. **Environment Agent**: Creates visualizations (dotplots) to examine marker gene expression patterns
4. **Evaluation Agent**: Analyzes visualizations and proposes cell-type annotations

This iterative process continues until satisfactory annotations are achieved.

## 📊 Input Data Requirements

Your AnnData object should contain:
- **Raw or normalized gene expression matrix** (`.X`)
- **Cell clustering results** (e.g., `.obs['leiden']` or `.obs['louvain']`)
- **Gene names** in `.var_names`

Optional but recommended:
- **UMAP coordinates** (`.obsm['X_umap']`)
- **PCA results** (`.obsm['X_pca']`)

## 🎯 Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_adata` | AnnData | None | Input AnnData object |
| `initial_hypothesis` | str | "" | Research context and expected cell types |
| `output_column` | str | "new_annotation" | Column name for annotations |
| `API_key` | str | "" | API key for LLM provider |
| `model_name` | str | "o1-mini" | LLM model to use |
| `model_provider` | str | "openai" | LLM provider |
| `original_grouping` | str | "leiden" | Clustering column name |
| `species` | str | "human" | Species (human/mouse) |
| `input_dir` | str | "CellMaster" | Input directory |
| `output_dir` | str | "CellMaster" | Output directory |

## 📈 Output

CellMaster generates:
- **Annotated AnnData object** with cell-type labels
- **UMAP plots** showing cell-type annotations
- **Dotplots** of marker gene expression
- **Log files** documenting the annotation process
- **Annotation CSV files** for each iteration

## 🧪 Example: PBMC Annotation

See `testing.ipynb` for a complete example. Here's a snippet:

```python
import scanpy as sc
from cli import run_cell_annotation

# Load PBMC data
adata = sc.read_h5ad("uploads/pbmc3k.h5ad")

# Define hypothesis
pbmc_hypothesis = """
You are a bioinformatics expert in single-cell RNA sequencing analysis.
Dataset: This is a scRNA dataset of 2638 cells and 1838 genes. It is a PBMC cell collection.
Project Description: This dataset contains PBMC cell types including T cells, B cells, 
Monocytes, NK cells, and their subgroups like CD4 T and CD8 T. You must include these cell types.
"""

# Run annotation
annotated_adata = run_cell_annotation(
    input_adata=adata,
    initial_hypothesis=pbmc_hypothesis,
    output_column="cellmaster_annotation",
    API_key=OPENAI_API_KEY,
    model_name="gpt-4o"
)

# View results
print(annotated_adata.obs["cellmaster_annotation"].value_counts())
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🔗 Related Projects

- **[CellMaster-UI](https://github.com/AnonymousGym/CellMaster-UI)**: Graphical user interface for CellMaster

## 📧 Contact

For questions, issues, or suggestions, please open an issue on GitHub.

## 🙏 Acknowledgments

CellMaster leverages several excellent open-source tools:
- [Scanpy](https://scanpy.readthedocs.io/) for single-cell analysis
- [AnnData](https://anndata.readthedocs.io/) for data structure
- OpenAI for LLM capabilities

## 📚 Citation

If you use CellMaster in your research, please cite our work (citation details coming soon).

---

**Note**: Always review automated annotations manually. CellMaster is a tool to assist researchers, not replace expert biological knowledge.