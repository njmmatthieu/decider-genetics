from enum import Enum, auto
import hashlib
import numpy as np
import pandas as pd
import random

class AllVariantsAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    PATIENT = auto()
    SAMPLE = auto()
    VARIANT = auto()


class AllVariantsAdapterPatientField(Enum):
    """
    Define possible fields the adapter can provide for patients.
    """

    ID = "patient"
    SAMPLES = "samples"


class AllVariantsAdapterSampleField(Enum):
    """
    Define possible fields the adapter can provide for samples.
    """

    ID = "sample"
    READ_COUNTS = "readCounts"


class AllVariantsAdapterVariantField(Enum):
    """
    Define possible fields the adapter can provide for variants.
    """

    ID = "ID"
    CHROMOSOME = "CHROM"
    POSITION = "POS"
    REF = "REF"
    ALT = "ALT"
    FILTER = "FILTER"
    CYTOBAND = "cytoBand"
    FUNCTION = "Func.MANE"
    GENE = "Gene.MANE"
    GENE_DETAIL = "GeneDetail.MANE"
    EXONIC_FUNCTION = "ExonicFunc.MANE"
    AA_CHANGE = "AAChange.MANE"
    FUNCTION_REF = "Func.refGene"
    GENE_REF = "Gene.refGene"
    GENE_DETAIL_REF = "GeneDetail.refGene"
    EXONIC_FUNCTION_REF = "ExonicFunc.refGene"
    AA_CHANGE_REF = "AAChange.refGene"
    GENOMIC_SUPER_DUPS = "genomicSuperDups"
    DBSC_SNV_ADA_SCORE = "dbscSNV_ADA_SCORE"
    DBSC_SNV_RF_SCORE = "dbscSNV_RF_SCORE"
    COSMIC_ID = "COSMIC_ID"
    COSMIC_OCCURENCE = "COSMIC_OCCURRENCE"
    COSMIC_TOTAL_OCCURENCE = "COSMIC_TOTAL_OCC"
    COSMIC_CONF_SOMATIC = "COSMIC_CONF_SOMA"
    CLNSIG = "CLNSIG"
    CLNSIGCONF = "CLNSIGCONF"
    CLNDN = "CLNDN"
    CLNREVSTAT = "CLNREVSTAT"
    CLNALLELEID = "CLNALLELEID"
    CLNDISDB = "CLNDISDB"
    INTERPRO_DOMAIN = "Interpro_domain"
    REGULOME_DB = "regulomeDB"
    CADD_RAW = "CADD_raw"
    CADD_PHRED = "CADD_phred"
    THOUSAND_GENOMES_ALL = "1000G_ALL"
    THOUSAND_GENOMES_EUR = "1000G_EUR"
    GNOMAD_GENOME_ALL = "gnomAD_genome_ALL"
    GNOMAD_GENOME_NFE = "gnomAD_genome_NFE"
    GNOMAD_GENOME_FIN = "gnomAD_genome_FIN"
    GNOMAD_GENOME_MAX = "gnomAD_genome_max"
    GNOMAD_EXOME_NC_ALL = "gnomAD_exome_nc_ALL"
    GNOMAD_EXOME_NC_NFE = "gnomAD_exome_nc_NFE"
    GNOMAD_EXOME_NC_NFE_SWE = "gnomAD_exome_nc_NFE_SWE"
    GNOMAD_EXOME_NC_NC_FIN = "gnomAD_exome_nc_FIN"
    GNOMAD_EXOME_NC_MAX = "gnomAD_exome_nc_max"
    TRUNCAL = "Truncal"


class AllVariantsAdapterEdgeType(Enum):
    """
    Enum for the edge types of the adapter.
    """

    PATIENT_SAMPLE_ASSOCIATION = auto()
    SAMPLE_VARIANT_ASSOCIATION = auto()
    VARIANT_GENE_ASSOCIATION = auto()

def preprocess_variants(df):

    df = pd.read_table("/Users/mnajm/Documents/DECIDER/Code/decider-genetics/data/synthetic_variants.csv", sep="\t")

    # break up the 'samples' column into one row per sample, rename the
    # column to 'sample'; at the same time, the readCounts and Gene.MANE
    # columns need to be split into one row per sample as well
    df = df.assign(
        samples=df["samples"].str.split(";")
    ).explode("samples")

    df = df.rename(
        columns={"samples": "sample"}
    ).reset_index(drop=True)

    df = df.assign(
        readCounts=df["readCounts"].str.split(";")
    ).explode("readCounts")

    df = df.assign(
        Gene=df["Gene.MANE"].str.split(";")
    ).explode("Gene")

    # remove duplicate rows
    df = df.drop_duplicates()

    # define dropped columns that should not contribute to hash
    _drop_columns = [
        AllVariantsAdapterPatientField.ID.value,
        AllVariantsAdapterSampleField.ID.value,
        AllVariantsAdapterSampleField.READ_COUNTS.value,
        "Gene",
    ]

    # if ID is '.', generate md5 hash from other columns
    df["ID"] = df.apply(
        lambda row: (
            hashlib.md5(
                "".join(
                    [
                        str(row[column])
                        for column in df.columns
                        if column not in _drop_columns
                    ]
                ).encode("utf-8")
            ).hexdigest()
            if row["ID"] == "."
            else row["ID"]
        ),
        axis=1,
    )

    return(df)
