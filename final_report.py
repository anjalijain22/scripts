import pandas as pd
import logging
import numpy as np


def log_message(*message):
    """write message to logfile and stdout"""
    if message:
        for i in message:
            logging.info(i)
            print(i)


def get_INFO_annot(df):
    """Return a dataframe with ClinVar and GnomAD annotations"""

    # dictionary to store the annotations we want to fetch from the INFO column
    data = {
        "clinvar_pathogenic": [],
        "gnomAD_AC_hom": [],
        "gnomAD_AC_het": [],
        "gnomAD_AF_hom": [],
        "gnomAD_AF_het": [],
        "gnomAD_max_hl": [],
    }

    for index, row in df.iterrows():
        INFO = row["INFO"].split(";")  # Access the INFO column and split by ";"
        INFO_names = [x.split("=")[0] for x in INFO]  # store the name of the info field
        INFO_values = [
            i.split("=")[1] for i in INFO
        ]  # store the value of the info field

        # If the key of data dict is present in INFO_names, then append corresponding INFO_value to the key's value,
        # else append value NA
        for key, value in data.items():
            if key in INFO_names:
                value.append(INFO_values[INFO_names.index(key)])
            else:
                value.append(".")

    # Create a df from the data dictionary
    annotations_df = pd.DataFrame(data)
    log_message(
        "Successfully fetched relevant ClinVar and GnomAD annotations from the INFO field of the report."
    )

    return annotations_df


def correct_mitotip_interpretations(df):
    """Change the MitoTip_interpretation column to correct values based on values of MitoTip_Score"""

    # Change datatype of column
    df["MitoTip_score"] = pd.to_numeric(df["MitoTip_score"], errors="coerce")

    df.loc[df["MitoTip_score"] > 8.44, "MitoTip_interpretation"] = "possibly benign"
    df.loc[
        df["MitoTip_score"] > 12.66, "MitoTip_interpretation"
    ] = "possibly pathogenic"
    df.loc[df["MitoTip_score"] > 16.25, "MitoTip_interpretation"] = "likely pathogenic"
    df["MitoTip_score"] = df["MitoTip_score"].replace(np.nan, ".")
    log_message("Corrected the MitoTip_interpretation column in report.")
    return df


def change_annot_9155(df):
    """Change status_mitomap column for m.9155A>G from . to Confirmed"""
    variant_entry = df.loc[df.HGVS == "m.9155A>G"]
    if variant_entry.empty:
        log_message("Variant m.9155A>G not present in the report.")
        return df
    else:
        df.loc[df.HGVS == "m.9155A>G", "status_mitomap"] = "Confirmed"
        log_message(
            "Found variant m.9155A>G in the report and updated status_mitomap to confirmed."
        )
        return df


def concat_df(df1, df2):
    """Concatenate two dataframes along axis 1 (column)"""

    concatenated_df = pd.concat([df1, df2], axis=1)
    log_message("Successfully joined the two dataframes")
    return concatenated_df


def remove_cols(df):
    """Remove unwanted columns from the report dataframe"""

    # List of column to be removed from the file
    remove_cols = [
        "total_sample_depth",
        "variant_quality",
        "locus_mitomap",
        "QUAL",
        "FILTER",
        "MQM_INFO",
        "MQMR_INFO",
        "QA_INFO",
        "QR_INFO",
        "SAF_INFO",
        "SAR_INFO",
        "SRF_INFO",
        "SRR_INFO",
        "SBR_INFO",
        "SBA_INFO",
        "POS_FILTER",
        "SBR_FILTER",
        "SBA_FILTER",
        "MQMR_FILTER",
        "AQR_FILTER",
        "GT_FORMAT",
        "QR_FORMAT",
        "AQR_FORMAT",
        "QA_FORMAT",
        "AQA_FORMAT",
        "INFO",
        "FORMAT",
    ]

    # remove columns and store the remaining cols in new_df
    new_df = df.drop(remove_cols, axis=1)
    log_message("Removed unwanted columns from the dataframe")
    return new_df


def sort_by_sample(df):
    subset_df = df[
        [
            "HGVS",
            "SAMPLE",
            "alt_depth",
            "ref_depth",
            "total_locus_depth",
            "variant_heteroplasmy",
        ]
    ]
    grouped_df = subset_df.groupby("HGVS").agg(
        {
            "SAMPLE": lambda x: ";".join(str(i) for i in x),
            "alt_depth": lambda x: ";".join(str(i) for i in x),
            "ref_depth": lambda x: ";".join(str(i) for i in x),
            "total_locus_depth": lambda x: ";".join(str(i) for i in x),
            "variant_heteroplasmy": lambda x: ";".join(str(i) for i in x),
        }
    )

    df = df.drop(
        [
            "SAMPLE",
            "alt_depth",
            "ref_depth",
            "total_locus_depth",
            "variant_heteroplasmy",
        ],
        1,
    )
    final = pd.merge(df, grouped_df, on="HGVS", how="outer")

    log_message("Report sorted by samples")
    return final.drop_duplicates(ignore_index=True)


def check_sort(df):
    sample = df.SAMPLE.unique()
    if len(sample) == 1:
        log_message("Only one sample present in report")
        return df
    else:
        log_message("Multiple samples present in report")
        updated_df = sort_by_sample(df)
        return updated_df


def reorder_cols(df):
    """Reorder columns in the report dataframe"""

    col_list = [
        "CHR",
        "POS",
        "REF",
        "ALT",
        "SAMPLE",
        "HGVS",
        "gene/locus",
        "gene/locus description",
        "variant_heteroplasmy",
        "alt_depth",
        "ref_depth",
        "total_locus_depth",
        "COHORT_COUNT",
        "disease_mitomap",
        "status_mitomap",
        "disease_amino_acid_change_mitomap",
        "allele_frequency_mitomap",
        "GenBank_frequency_mitomap",
        "clinvar_pathogenic",
        "gnomAD_AC_hom",
        "gnomAD_AC_het",
        "gnomAD_AF_hom",
        "gnomAD_AF_het",
        "gnomAD_max_hl",
        "homoplasmy_mitomap",
        "heteroplasmy_mitomap",
        "num_references_mitomap",
        "variant_amino_acid_change_mitomap",
        "codon_position_mitomap",
        "codon_number_mitomap",
        "num_disease_references_mitomap",
        "RNA_mitomap",
        "tier",
        "commercial_panels",
        "phylotree_haplotype",
        "MitoTip_score",
        "MitoTip_interpretation",
        "MitoTip_percentile",
        "anticodon",
        "MGRB_frequency",
        "MGRB_FILTER",
        "MGRB_AC",
        "MGRB_AN",
        "phylotree_mut",
    ]

    reordered_df = df[col_list]
    log_message("Rearanged the columns in the dataframe")
    return reordered_df


def main(report):
    # family = snakemake.wildcards.family
    # logfile =  f"logs/mity/report/{family}_final_report_generation.log"
    logfile = f"/hpf/largeprojects/ccmbio/ajain/mity/test2_mity/test_out/ashkenazim.report_gen_2.log"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )

    report_df = pd.read_csv(report)

    CV_gnomad_annots_df = get_INFO_annot(report_df)
    merged_report = concat_df(report_df, CV_gnomad_annots_df)
    final_report = remove_cols(merged_report)
    final_report = check_sort(final_report)
    final_report = reorder_cols(final_report)
    final_report = correct_mitotip_interpretations(final_report)
    final_report = change_annot_9155(final_report)

    # final_report.to_csv("{family}_mity_vcfanno_final_report.csv", index=False)
    final_report.to_csv(
        "/hpf/largeprojects/ccmbio/ajain/mity/test2_mity/test_out/ashkenazim.vcfanno.processed.Apr18_2.csv",
        index=False,
    )

    log_message(
        "Final formatted report containing annotated list of mitochondrial variants created!"
    )


if __name__ == "__main__":
    report = "/hpf/largeprojects/ccmbio/ajain/mity/test2_mity/test_out/ashkenazim.annotated_variants.csv"
    # report = snakemake.input.report
    main(report)
