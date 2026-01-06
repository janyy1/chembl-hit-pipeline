import pandas as pd
import numpy as np


# cleaning
def remove_non_numeric_activity(df: pd.DataFrame) -> pd.DataFrame:
    # 숫자가 아닌 standard_value 제거
    
    df_clean = df.copy()
    
    df_clean['standard_value'] = pd.to_numeric(df_clean['standard_value'], errors='coerce')

    df_clean = df_clean.dropna(subset=['standard_value'])

    return df_clean


def convert_activity_to_nM(df: pd.DataFrame) -> pd.DataFrame:
    # standard_value 단위 nM으로 통일
    
    df_convert = df.copy()

    UNIT_TO_NM = {'nM': 1, 'uM': 1e3, 'µM': 1e3, 'pM': 1e-3, 'M': 1e9}

    df_convert['unit_factor'] = df_convert['standard_units'].map(UNIT_TO_NM) # 배율에 해당하는 임시 컬럼 만들기
    df_convert = df_convert.dropna(subset=['unit_factor']) # 변환 불가능한 값 제거
    
    df_convert['standard_value'] = (df_convert['standard_value'] * df_convert['unit_factor']) # standard_value 값에 배율을 곱함
    df_convert['standard_units'] = 'nM' # 단위 전체 nM으로 변환
    
    df_convert = df_convert.drop(columns=['unit_factor'])

    return df_convert


def remove_invalid_activity_range(df: pd.DataFrame) -> pd.DataFrame:
    # 숫자이긴 하지만, 생물학적으로 의미 없거나 측정 오류가 의심되는 활성값 범위 제거

    df_remove = df.copy()

    lower_bound = 0.1
    upper_bound = 10000

    condition = (
        (df_remove['standard_value'] > lower_bound) &
        (df_remove['standard_value'] < upper_bound)
        )

    df_remove = df_remove.loc[condition]

    return df_remove


# filtering
def filter_confidence_score(df: pd.DataFrame, min_score: int = 7) -> pd.DataFrame:
    # 신뢰도 낮은 activity 제거

    if 'confidence_score' not in df.columns:
        return df

    return df[df['confidence_score'] >= min_score]

    
def filter_assay_type(df: pd.DataFrame) -> pd.DataFrame:
    # ADME 거르기

    return df[df['assay_type'].isin(['B', 'F'])]


def filter_IC50(df: pd.DataFrame) -> pd.DataFrame:
    # (필요한 경우) IC50 데이터만 남기기

    return df[df['standard_type'] == 'IC50']



# aggregation
def aggregate_duplicate_activity(df: pd.DataFrame, method: str = 'median') -> pd.DataFrame:
    # 동일 화합물의 multiple measurement 통합

    df_aggre = df.copy()
    std = ['molecule_chembl_id', 'standard_type']

    df_aggre = df_aggre.groupby(std, as_index=False).agg(standard_value=('standard_value', 'median'), pIC50=('pIC50', 'median'), pKi=('pKi', 'median'))

    return df_aggre


# transformation
def add_pActivity(df: pd.DataFrame) -> pd.DataFrame:
    # 활성을 pX 값으로 변환

    df = df.copy()

    positive_mask = df['standard_value'] > 0
    df = df.loc[positive_mask]

    df['pActivity'] = -np.log10(df['standard_value'] * 1e-9) # pActivity 식 계산

    df['pIC50'] = np.where(
        df['standard_type'] == 'IC50',
        df['pActivity'],
        np.nan
        )

    df['pKi'] = np.where(
        df['standard_type'] == 'Ki',
        df['pActivity'],
        np.nan
        )    
    
    return df


# processed df 최종 완성
def preprocess_activity_data(df: pd.DataFrame) -> pd.DataFrame:
    df = remove_non_numeric_activity(df)
    df = convert_activity_to_nM(df)
    df = remove_invalid_activity_range(df)
    df = add_pActivity(df)
    df = filter_IC50(df)
    df = filter_assay_type(df)
    df = filter_confidence_score(df)
    print("Before aggregation:", df.columns.tolist())
    df = aggregate_duplicate_activity(df)
    print("After aggregation:", df.columns.tolist())

    required_cols = [
        'molecule_chembl_id',
        'canonical_smiles',
        'standard_type',
        'standard_value',
        'standard_units',
        'assay_type',
        'confidence_score',
        'pActivity',
        'pIC50',
        'pKi'
        ]

    df = df[[c for c in required_cols if c in df.columns]]

    return df


