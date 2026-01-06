import pandas as pd
from .chembl_loader import ChEMBLLoader
from .preprocessor import preprocess_activity_data

class Analyzer(ChEMBLLoader):

    def __init__(self, df):
        self.df = preprocess_activity_data(df)

    def identify_hits(
        self,
        cutoff: float = 6.0,
        activity_col: str = 'pIC50',
        assay_type: str | None = 'B',
        min_n: int = 1
        ):

        df_hits = self.df.copy()

        df_hits = self.df[self.df['pIC50'].notna()]
        df_hits = df_hits[df_hits['pIC50'] >= cutoff]

        counts = df_hits.groupby('molecule_chembl_id').size() # measurement 수 계산
        valid_compounds = counts[counts >= min_n].index # counts가 min_n 을 넘기는 compound만 남겨서, id 목록으로
        df_hits = df_hits[df_hits['molecule_chembl_id'].isin(valid_compounds)] # id 목록에 포함된 행만 유지

        if df_hits.empty:
            return df_hits
        
        #SettingWithCopyWarning
        df_hits = df_hits.copy()
        df_hits['is_hit'] = True

        return df_hits


    # identify_hits()가 반환한 row-level hit 데이터를 compound-level 의사결정용 테이블로 바꾸는 함수
    def summarize_hits(
        self,
        df_hits: pd.DataFrame,
        activity_col: str = 'pIC50'
        ) -> pd.DataFrame:

        if df_hits.empty:
            return pd.DataFrame(
                columns=[
                    'molecule_chembl_id',
                    'measurement_count',
                    'median_activity',
                    'mean_activity',
                    'std_activity',
                    'min_activity',
                    'max_activity'
                    ]
                )

        summary = (
            df_hits
            .groupby('molecule_chembl_id')
            .agg(
                measurement_count = (activity_col, 'count'),
                median_activity = (activity_col, 'median'),
                mean_activity = (activity_col, 'mean'),
                std_activity = (activity_col, 'std'),
                min_activity = (activity_col, 'min'),
                max_activity = (activity_col, 'max')
                )
            .reset_index()
            )

        return summary


    def classify_hit_strength(
        self,
        df_summary: pd.DataFrame,
        strong_cutoff: float = 7.0,
        weak_cutoff: float = 6.0,
        max_std: float = 1.5,
        min_n: int = 1
        ) -> pd.DataFrame:

        df = df_summary.copy()

        df['hit_strength'] = 'non-hit'

        df['pass_activity'] = df['median_activity'] >= strong_cutoff
        df['pass_std'] = df['std_activity'].notna() & (df['std_activity'] <= max_std)
        df['pass_n'] = df['measurement_count'] >= min_n

        strong_mask = df['pass_activity'] & df['pass_std'] & df['pass_n']

        weak_mask = (
            (df['median_activity'] >= weak_cutoff) &
            (~strong_mask)
            )

        df.loc[strong_mask, 'hit_strength'] = 'strong'
        df.loc[weak_mask, 'hit_strength'] = 'weak'

        ambiguous_mask = (
            (df['median_activity'] >= weak_cutoff) &
            (df['hit_strength'] == 'non-hit')
            )

        df.loc[ambiguous_mask, 'hit_strength'] = 'ambiguous'

        return df
