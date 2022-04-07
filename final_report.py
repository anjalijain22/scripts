import pandas as pd
import logging


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
        "gnomad_mitotip_score": [],
        "gnomad_mitotip_trna_prediction": [],
        "AC_hom": [],
        "AC_het": [],
        "AF_hom": [],
        "AF_het": [],
        "max_hl": [],
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
        "MitoTip_score",
        "MitoTip_interpretation",
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
        "tier",
        "commercial_panels",
        "phylotree_haplotype",
        "gnomad_mitotip_score",
        "gnomad_mitotip_trna_prediction",
        "MitoTip_percentile",
        "anticodon",
        "allele_frequency_mitomap",
        "disease_mitomap",
        "MGRB_frequency",
        "MGRB_FILTER",
        "MGRB_AC",
        "MGRB_AN",
        "phylotree_mut",
        "locus_mitomap",
        "num_references_mitomap",
        "variant_amino_acid_change_mitomap",
        "codon_position_mitomap",
        "codon_number_mitomap",
        "num_disease_references_mitomap",
        "RNA_mitomap",
        "homoplasmy_mitomap",
        "heteroplasmy_mitomap",
        "status_mitomap",
        "disease_amino_acid_change_mitomap",
        "GenBank_frequency_mitomap",
        "clinvar_pathogenic",
        "AC_hom",
        "AC_het",
        "AF_hom",
        "AF_het",
        "max_hl",
    ]

    reordered_df = df[col_list]
    log_message("Rearanged the columns in the dataframe")
    return reordered_df


def main(report):
    family = snakemake.wildcards.family
    logfile =  f"logs/mity/report/{family}_final_report_generation.log"
    #logfile = f"1760_rep.log"
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
    final_report = reorder_cols(final_report)

    final_report.to_csv("{family}_mity_vcfanno_final_report.csv", index=False)
    #final_report.to_csv("1760_mity_vcfanno_final_report2.csv", index=False)

    log_message(
        "Final formatted report containing annotated list of mitochondrial variants created!"
    )


if __name__ == "__main__":
    #report = "1760_mito.annotated_variants.csv"
    report = snakemake.input.report
    main(report)