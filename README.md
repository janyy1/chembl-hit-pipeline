# ChEMBL Hit Identification Pipeline

This project is a **Python-based data analysis pipeline** designed to process public bioactivity data from the ChEMBL database and systematically identify potential hit compounds from noisy and duplicated experimental measurements.

The main goals of this project are:

- Practicing data processing using **pandas**

- Designing an analysis pipeline with an **OOP-based structure**

- Understanding the practical data flow of **drug discovery and bioactivity analysis**

This is an **educational project**, focusing on analysis logic and decision-making rather than producing validated biological conclusions.

---

## Project Goals

- Understand the structure of ChEMBL bioactivity data

- Aggregate duplicated experimental results into reliable representative values

- Define hit compounds based on clear and adjustable rules

- Design analysis logic as **reusable and maintainable Python classes**

---

## Project Scope

This project focuses on building a reproducible data analysis pipeline for hit identification 
based on public bioactivity data. It does not aim to develop predictive machine learning models 
or to replace experimental validation.

---

## Project Structure

```
chembl-hit-pipeline/
│
├── README.md
├── requirements.txt
│
├── src/
│   └── chembl_pipeline/
│       ├── __init__.py
│       ├── loader.py        # ChEMBLLoader: data loading
│       ├── preprocess.py   # preprocessing functions
│       ├── analyzer.py     # Analyzer: hit identification and summary
│
├── notebooks/
│   └── exploration.ipynb   # exploratory data analysis
│
└── output/
│   └── results.json        # example analysis output
```

---

## Data Overview

- Source: ChEMBL

- Main columns used
    - molecule_chembl_id
    - canonical_smiles
    - standard_type (e.g. IC50)
    - standard_value
    - standard_units
    - assay_chembl_id

This project focuses on **quantitative bioactivity values** (such as IC50) for downstream analysis.

---

## Workflow

### 1. Data Loading

- Load ChEMBL activity data into a pandas DataFrame

- Select only the columns required for analysis

### 2. Preprocessing

- Remove missing or invalid values

- Filter by activity type and units

- Aggregate duplicated activity measurements (same compound under the same conditions)

By default, **median aggregation** is used to reduce the influence of outliers,
though other methods (e.g. mean) can be selected depending on the analysis goal.

### 3. Hit Identification

- Identify hit candidates based on predefined rules

- Example criteria:
    - IC50 ≤ a specified threshold
    - Minimum number of experimental measurements

The threshold values are based on ranges commonly used in the literature and
are implemented as **configurable parameters** to allow flexible adjustment.

### 4. Result Summary

- Generate representative activity values for each compound

- Classify compounds as hit / non-hit

- Save results in JSON format

---

## Design Notes

- Object-Oriented Design
    - ChEMBLLoader: responsible for data loading
    - Analyzer: responsible for analysis logic and hit criteria

This separation allows changes to data sources or hit definitions without modifying the entire pipeline.

- Preprocessing logic is implemented as standalone functions to improve:
    - Testability
    - Reusability
    - Flexibility when analysis criteria change

---

## Example Usage

```markdown
```python
from chembl_pipeline.loader import ChEMBLLoader
from chembl_pipeline.analyzer import Analyzer

loader = ChEMBLLoader("data/chembl_activity.csv")
df = loader.load()

analyzer = Analyzer(df, hit_threshold=1000)
hit_results = analyzer.identify_hits()

analyzer.save_json("output/results.json")
```

---

## Limitations

This pipeline aggregates multiple experimental measurements into representative values
to enable compound-level analysis. As a result, differences in individual assay conditions
(e.g., experimental setup or protocol variability) are not fully reflected.

The identified hits should be interpreted as **screening-level candidates** for prioritization,
not as validated leads or absolute activity measurements.