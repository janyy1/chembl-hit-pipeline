import pandas as pd
from chembl_webresource_client.new_client import new_client

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)

class ChEMBLLoader:

    def __init__(
        self,
        target_chembl_id: str,
        assay_type=None,
        debug=False
        ):
        self.target_chembl_id = target_chembl_id
        self.assay_type = assay_type
        self.debug = debug
        self.raw_data = None
        self.df = None


    def fetch_bioactivities(self, standard_types=None):

        activity = new_client.activity

        query = activity.filter(
            target_chembl_id=self.target_chembl_id)

        if standard_types:
            query = query.filter(standard_type__in=standard_types)

        if self.debug:
            query = query.filter(confidence_score__gtes=8)
            query = query.filter(assay_type='B')
            
        self.raw_data = list(query.only([
                'molecule_chembl_id',
                'canonical_smiles',
                'standard_type',
                'standard_value',
                'standard_units',
                'assay_type',
                'confidence_score'
                ]))[:200]
        
        return self.raw_data


    def to_dataframe(self):
        if self.raw_data is None:
            raise ValueError("fetch_bioactivities() 먼저 실행하세요")

        self.df = pd.DataFrame(list(self.raw_data))
        return self.df
