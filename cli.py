import os
from pathlib import Path
from typing import Dict, Any
from pipeline import CellAnnotationPipeline

def run_cell_annotation(
    *,
    model_name: str = "o1-mini",
    model_provider: str = "openai",
    API_key = "",
    input_dir: str = "CellMaster",
    output_dir: str = "CellMaster",
    log_file_path: str = "CellMaster/annotation_log_file.txt",
    h5ad_file: str = "CellMaster.h5ad",
    markers_file: str = "CellMaster.csv",
    original_grouping: str = "leiden",
    output_column: str = "new_annotation",
    species: str = "human",
    initial_hypothesis: str = "",
    input_adata= None,
    **kwargs,
):
    """
    Run the cell-annotation pipeline and return the annotated AnnData object.
    """     
    # 1. Normalise / validate paths
    input_dir      = Path(input_dir) # type: ignore
    output_dir     = Path(output_dir) # type: ignore
    log_file_path  = Path(log_file_path) # type: ignore
    h5ad_file      = Path(h5ad_file) # type: ignore

    output_dir.mkdir(parents=True, exist_ok=True) # type: ignore
    log_file_path.parent.mkdir(parents=True, exist_ok=True) # type: ignore

    if input_adata:
        input_adata.write(os.path.join(input_dir, h5ad_file))

    # 2. Assemble the config expected by `CellAnnotationPipeline`
    config: Dict[str, Any] = {
        "model_name":        model_name,
        "model_provider":    model_provider,
        "API_key":           API_key,
        "input_dir":         str(input_dir),
        "output_dir":        str(output_dir),
        "log_file_path":     str(log_file_path),
        "h5ad_file":         str(h5ad_file),
        "markers_file":      str(markers_file),
        "original_grouping": original_grouping,
        "output_column":     output_column,
        "species":           species,
        "initial_hypothesis": initial_hypothesis,
        # Allow callers to override / extend:
        **kwargs,
    }

    ITERATION = 1
    MAX_RETRIES = 3

    # 3. Instantiate and run the pipeline
    pipeline = CellAnnotationPipeline(config)
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            pipeline.run_pipeline(iterations=ITERATION)
            break  # Success, so break the loop
        except Exception as e:
            print(f"Attempt {attempt} failed with error: {e}")
            if attempt < MAX_RETRIES:
                print("Retrying...")
            else:
                print("All retries failed.")
                
    # 4. Return the annotated AnnData
    return pipeline.adata
