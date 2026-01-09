
import os
import logging
from typing import Optional, List, Literal
import pandas as pd
from pathlib import Path

from matrix_eqtl.preprocessing import genotype, expression, covariates
from matrix_eqtl.tools import matrix_eqtl, qtltools, tensorqtl
from matrix_eqtl.utils import io

class EQTLPipeline:
    def __init__(self, output_dir: str, tool: Literal['matrixeqtl', 'qtltools', 'tensorqtl'] = 'matrixeqtl'):
        """
        Initialize the eQTL analysis pipeline.
        
        Args:
            output_dir: Directory to store output files.
            tool: The analysis tool to use ('matrixeqtl', 'qtltools', 'tensorqtl').
        """
        self.output_dir = Path(output_dir)
        self.tool_name = tool
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        logging.basicConfig(level=logging.INFO, 
                            format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)

    def run(self, 
            vcf_file: str, 
            expression_file: str, 
            covariates_file: Optional[str] = None,
            gene_location_file: Optional[str] = None,
            snp_location_file: Optional[str] = None,
            genotype_pca_n: int = 3,
            peer_factors_n: int = 0,
            cis_window: int = 1000000,
            output_location_matrix: bool = False):
        """
        Execute the full pipeline.
        
        Args:
            vcf_file: Path to input VCF file.
            expression_file: Path to expression matrix.
            covariates_file: Path to additional covariates file.
            gene_location_file: Path to gene location file (needed for cis-eQTL).
            genotype_pca_n: Number of genotype PCs to calculate.
            peer_factors_n: Number of PEER factors to calculate.
            cis_window: Window size for cis-eQTL analysis.
            output_location_matrix: Whether to output formatted location matrices.
        """
        self.logger.info(f"Starting pipeline using tool: {self.tool_name}")
        
        # 0. Intersect Samples
        self.logger.info("Intersecting samples...")
        common_samples = io.get_common_samples(vcf_file, expression_file)
        self.logger.info(f"Fouund {len(common_samples)} common samples.")
        
        # 1. Process Expression (Filter -> Norm)
        self.logger.info("Processing expression data...")
        processed_exp = expression.process_expression(
            expression_file=expression_file,
            samples=common_samples,
            output_dir=self.output_dir / "preprocessing",
            filter_low_expr=True,
            normalize_iqn=True
        )
        
        # 2. Process Genotype
        self.logger.info("Processing genotype data...")
        processed_geno = genotype.process_genotype(
            vcf_file=vcf_file,
            samples=common_samples,
            output_dir=self.output_dir / "preprocessing",
            target_tool=self.tool_name
        )
        
        # 3. Generate/Merge Covariates
        self.logger.info("Generating covariates...")
        final_covariates = covariates.generate_covariates(
            genotype_prefix=processed_geno['plink_prefix'],
            expression_matrix=processed_exp,
            known_covariates=covariates_file,
            output_dir=self.output_dir / "covariates",
            pca_n=genotype_pca_n,
            peer_n=peer_factors_n
        )
        
        # 4. Run Analysis Tool
        self.logger.info(f"Running {self.tool_name}...")
        
        if self.tool_name == 'matrixeqtl':
            matrix_eqtl.run_analysis(
                genotype_file=processed_geno['matrix'],
                expression_file=processed_exp,
                covariates_file=final_covariates,
                output_dir=self.output_dir / "results",
                gene_loc=gene_location_file,
                snp_loc=processed_geno.get('snp_loc'),
                cis_window=cis_window
            )
        elif self.tool_name == 'qtltools':
            qtltools.run_analysis(
                vcf_file=processed_geno['vcf'], # QTLtools uses VCF/BED
                expression_bed=processed_exp,    # Needs BED format
                covariates_file=final_covariates,
                output_dir=self.output_dir / "results",
                gene_loc=gene_location_file, # Embedded in BED usually
                cis_window=cis_window
            )
        elif self.tool_name == 'tensorqtl':
            tensorqtl.run_analysis(
                plink_prefix=processed_geno['plink_prefix'],
                expression_file=processed_exp,
                covariates_file=final_covariates,
                output_dir=self.output_dir / "results",
                cis_window=cis_window
            )
            
        self.logger.info("Pipeline completed successfully.")
