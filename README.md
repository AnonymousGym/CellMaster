# CellMaster

**CellMaster** is an intelligent automated pipeline for single-cell RNA sequencing (scRNA-seq) cell-type annotation. It leverages Large Language Models (LLMs) to iteratively refine cell-type annotations through a multi-agent system that includes hypothesis generation, marker gene proposal, visualization, and evaluation.

## Features

- **Automated Cell-Type Annotation**: Intelligent, iterative annotation of cell types in scRNA-seq data
- **LLM-Powered Analysis**: Utilizes advanced language models (gpt-4o, o1) for biological reasoning
- **Multi-Agent Architecture**: Coordinated agents for hypothesis, experiment, environment, and evaluation
- **Visualization**: Automatic generation of UMAP plots and dotplots for validation
- **Iterative Refinement**: Progressive improvement of annotations through multiple iterations
- **Marker Gene Discovery**: Intelligent proposal and validation of cell-type-specific marker genes

## Quick Start

### Demo / Vignette

**See [`testing.ipynb`](testing.ipynb) for a complete demonstration of CellMaster's capabilities.** This Jupyter notebook provides a step-by-step example of cell-type annotation on demo datasets.

### UI Version

For a graphical user interface version of CellMaster, check out **[@AnonymousGym/CellMaster-UI](https://github.com/AnonymousGym/CellMaster-UI)**.

## Prerequisites

- Python 3.11+
- OpenAI API key or compatible LLM API access
- Single-cell RNA-seq data in `.h5ad` format (AnnData)

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/AnonymousGym/CellMaster.git
cd CellMaster
```

### 2. Set Up Environment

Install dependencies:

```bash
pip install -r requirements.txt
```

### 3. Configure API Key

**Recommended: Use Environment Variables (Most Secure)**

Set your API key as an environment variable:

```bash
export OPENAI_API_KEY='your-api-key-here'
```

Or create a `.env` file in the project root (and add `.env` to `.gitignore`):

```bash
# .env
OPENAI_API_KEY=your-api-key-here
```

Then load it in Python:

```python
import os
from dotenv import load_dotenv  # pip install python-dotenv

load_dotenv()  # Load environment variables from .env file
api_key = os.environ.get('OPENAI_API_KEY')
```

### 4. Download Required Data Files

Download the example dataset package from Google Drive:
- [Example Datasets (retina.h5ad, pbmc3k.h5ad)](https://drive.google.com/file/d/1qT5qBbUYP3BvXNL7NY6zqEr-pp13zsdW/view?usp=sharing)

Place the downloaded files in an `uploads/` directory in the project root.

## Usage

### Basic Usage with Python API

```python
from cli import run_cell_annotation
import scanpy as sc
import os

# Load your scRNA-seq data
adata = sc.read_h5ad("path/to/your/data.h5ad")

# Get API key from environment variable (recommended)
api_key = os.environ.get('OPENAI_API_KEY')

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
    API_key=api_key,
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
    API_key=os.environ.get('OPENAI_API_KEY'),  # Use environment variable
    model_name="gpt-4o",
    model_provider="openai",
    original_grouping="leiden",  # clustering column name
    species="human",  # or "mouse"
    input_dir="CellMaster",
    output_dir="CellMaster",
    log_file_path="CellMaster/annotation_log.txt"
)
```

## Input Data Requirements

Your AnnData object should contain:
- **Raw or normalized gene expression matrix** (`.X`)
- **Cell clustering results** (e.g., `.obs['leiden']` or `.obs['louvain']`)
- **Gene names** in `.var_names`

## Parameters

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

## Output

CellMaster generates:
- **Annotated AnnData object** with cell-type labels
- **UMAP plots** showing cell-type annotations
- **Dotplots** of marker gene expression
- **Log files** documenting the annotation process
- **Annotation CSV files** for each iteration

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Related Projects

- **[CellMaster-UI](https://github.com/AnonymousGym/CellMaster-UI)**: Graphical user interface for CellMaster

---

**Note**: Always review automated annotations manually. CellMaster is a tool to assist researchers, not replace expert biological knowledge.
