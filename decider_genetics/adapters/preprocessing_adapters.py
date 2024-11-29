import pandas as pd
import random

import ontoweaver


class Clinical(ontoweaver.tabular.PandasAdapter):

    def __init__(self,
        df: pd.DataFrame,
        config: dict,
    ):
        
        df = {k.lower().replace(" ", "_"): v for k, v in df.items()}
        df["bmi"] = float(df["bmi"].replace(",", "."))
        df["chemotherapy_cycles"] = int(df["chemotherapy_cycles"])

        # convert parpi, brca_mutation, and hr_deficient to boolean
        df["parpi"] = str(df["parpi"]).lower() == "yes"
        df["brca_mutation"] = (
            str(df["brca_mutation"]).lower() == "yes"
        )
        df["hr_deficient"] = (
            str(df["hr_deficient"]).lower() == "hrd positive"
        )
        # convert bool to lowercase string
        df["parpi"] = str(df["parpi"]).lower()
        df["brca_mutation"] = str(df["brca_mutation"]).lower()
        df["hr_deficient"] = str(df["hr_deficient"]).lower()

        # add fake severe_adverse_reaction randomly
        drugs = [
            "cisplatin",
            "bevacizumab",
            "abeciclib",
            "olaparib",
            "paclitaxel",
        ]
        # with 20% probability, sample one of the drugs
        if random.random() < 0.2:
            df["severe_adverse_reaction_to"] = random.choice(drugs)

        # Default mapping as a simple config.
        from . import types
        parser = ontoweaver.tabular.YamlParser(config, types)
        mapping = parser()

        # Declare types defined in the config.
        super().__init__(
            df,
            *mapping,
        )