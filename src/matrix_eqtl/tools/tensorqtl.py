
import subprocess
from pathlib import Path

def run_analysis(plink_prefix, expression_file, covariates_file, output_dir, cis_window=1e6):
    """
    Run TensorQTL analysis.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / "tensorqtl_run"
    
    # Assuming we wrap the vendor/tensorqtl or use the wrapper script I moved
    script_path = Path(__file__).parent.parent / "analysis" / "run_tensorqtl.py"

    cmd = [
        "python3", str(script_path),
        "--plink", str(plink_prefix),
        "--expression", str(expression_file),
        "--covariates", str(covariates_file),
        "--output", str(prefix),
        "--window", str(cis_window)
    ]
    
    print(f"Running TensorQTL: {' '.join(cmd)}")
    subprocess.check_call(cmd)
