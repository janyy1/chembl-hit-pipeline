from chembl_pipeline.chembl_loader import ChEMBLLoader
from chembl_pipeline.analyzer import Analyzer
from pathlib import Path

def main():

    # 1. Load data
    loader = ChEMBLLoader(
        target_chembl_id="CHEMBL204",
        assay_type="B"
        )

    print("Fetching bioactivities...")
    loader.fetch_bioactivities(standard_types=["IC50", "Ki"])
    print("Converting to DataFrame...")
    df_raw = loader.to_dataframe()
    print("DataFrame shape:", df_raw.shape)

    # 2. Analyze
    analyzer = Analyzer(df_raw)

    df_hits = analyzer.identify_hits(
        cutoff=6.0,
        activity_col="pIC50",
        assay_type="B",
        min_n=1
        )

    df_summary = analyzer.summarize_hits(df_hits)
    df_summary = analyzer.classify_hit_strength(df_summary)

    # 3. Output
    print("=== Hit Summary ===")
    print(df_summary.head())

    output_dir = Path("results")
    output_dir.mkdir(exist_ok=True)
    
    df_summary.to_csv(output_dir / "hit_summary.csv", index=False)
    


if __name__ == "__main__":
    main()
