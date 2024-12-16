from enum import Enum
import hashlib
import numpy as np
import pandas as pd
import random


class CnGenesAdapterSampleField(Enum):
    """
    Define possible fields the adapter can provide for samples.
    """

    ID = "sample"


class CnGenesAdapterGeneField(Enum):
    """
    Define possible fields the adapter can provide for genes.
    """

    ENSEMBL_ID = "ID"
    NAME = "Gene"
    CHR = "chr"
    START = "start"
    END = "end"
    STRAND = "strand"
    BAND = "band"
    TYPE = "type"  # should be prefiltered to only include 'protein_coding'

class CnGenesAdapterEdgeField(Enum):
    """
    Enum for the edge fields of the adapter.
    """

    N_PROBES_CR = "nProbesCr"
    N_PROBES_AF = "nProbesAf"
    LOG_R = "logR"
    BAF = "baf"
    N_ARAW = "nAraw"
    N_BRAW = "nBraw"
    N_MAJOR = "nMajor"
    N_MINOR = "nMinor"
    PURIFIED_LOG_R = "purifiedLogR"
    PURIFIED_BAF = "purifiedBaf"
    PURIFIED_LOH = "purifiedLoh"
    CN_STATUS = "CNstatus"
    LOH_STATUS = "LOHstatus"
    MIN_PURIFIED_LOG_R = "minPurifiedLogR"
    MAX_PURIFIED_LOG_R = "maxPurifiedLogR"
    BREAKS_IN_GENE = "breaksInGene"

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
        return np.nan

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

    df['severe_adverse_reaction_to']= [generate_row(drugs) for _ in range(df.shape[0])]

    return df

def df_to_string(t):
    # Convert Tuple to String
    res = '_'.join(str(val) for val in t)
    return res

def preprocess_cna(df):

    columns_to_id = [field.value for field in CnGenesAdapterEdgeField]+[CnGenesAdapterGeneField.NAME.value,CnGenesAdapterSampleField.ID.value]
    df['VARIANT_ID'] = df[columns_to_id].apply(tuple, axis=1).apply(df_to_string)

    # # Remove all columns except the ones in
    # # CnGenesAdapterEdgeField, plus the sample (patient) id and gene NAME columns, and
    # # deduplicate
    # df = df[
    # [
    #     field.value
    #     for field in CnGenesAdapterEdgeField
    #     if field.value in df.columns
    # ]
    # + [
    #     CnGenesAdapterSampleField.ID.value,
    #     CnGenesAdapterGeneField.NAME.value,
    #     ':ID'
    # ]
    # ].drop_duplicates()

    # # generate an id for each variant using the md5 hash of all columns
    # df["VARIANT_ID"] = df.apply(
    #     lambda row: hashlib.md5(
    #         "".join(
    #             [str(row[column]) for column in df.columns]
    #         ).encode("utf-8")
    #     ).hexdigest(),
    #     axis=1,
    # )

    return(df)