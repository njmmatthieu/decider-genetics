import pandas as pd
import random

# import ontoweaver

def preprocess_clinical(df):
        
    # lowercase the keys, replace space with underscore
    df.columns = map(str.lower, df.columns)
    df.columns = df.columns.str.replace(pat=" ", repl="_", regex=True)
    print(df.columns)

    df["bmi"] = df["bmi"].str.replace(pat=",", repl=".", regex=True).astype(float)
    df["chemotherapy_cycles"] = df["chemotherapy_cycles"].astype(int)

    # convert parpi, brca_mutation, and hr_deficient to boolean
    df["parpi"] = df["parpi"].str.lower() == "yes"
    df["brca_mutation"] = df["brca_mutation"].str.lower() == "yes"
    df["hr_deficient"] = df["hr_deficient"].str.lower() == "hrd positive"

    # convert bool to lowercase string
    df["parpi"] = df["parpi"].astype(str).str.lower()
    df["brca_mutation"] = df["brca_mutation"].astype(str).str.lower()
    df["hr_deficient"] = df["hr_deficient"].astype(str).str.lower()

    # add fake severe_adverse_reaction randomly
    drugs = [
        "cisplatin",
        "bevacizumab",
        "abeciclib",
        "olaparib",
        "paclitaxel",
    ]

    # add fake severe_adverse_reaction randomly
    def generate_row(drugs):
        """
        With 20% probability, sample one of the drugs.

        Args:
            drugs (list): The list of drugs.

        Returns:
            str: The severe adverse reaction to a randomly chosen drug.
        """
        if random.random() < 0.2:
            return random.choice(drugs)
        else:
            return None

    df['severe_adverse_reaction_to']= [generate_row(drugs) for _ in range(df.shape[0])]

    return df