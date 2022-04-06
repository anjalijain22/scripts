import pandas as pd

#file to be formatted
file="1760_mito.annotated_variants.csv"

#import file as pandas dataframe
file_df=pd.read_csv(file)

#dictionary to store the annotations we want to fetch from the INFO column
data={"clinvar_pathogenic":[], "gnomad_mitotip_score":[], "gnomad_mitotip_trna_prediction":[], 
"AC_hom":[], "AC_het":[], "AF_hom":[], "AF_het":[], "max_hl":[]}


for index,row in file_df.iterrows():
    INFO=row["INFO"].split(";") #Access the INFO column and split by ";"
    INFO_names = [x.split("=")[0] for x in INFO] #store the name of the info field
    INFO_values = [i.split("=")[1] for i in INFO] #store the value of the info field
    
    #If the key of data dict is present in INFO_names, then append corresponding INFO_value to key's value,
    # else append value NA
    for key,value in data.items():
        if key in INFO_names:
            value.append(INFO_values[INFO_names.index(key)])
        else:
            value.append(".")

#Create a df from the data dictionary
CV_gnomad_annot=pd.DataFrame(data)
print(CV_gnomad_annot)

#FINAL FORMATTING OF THE REPORT

#List of column to be removed from the file
remove_cols=["total_sample_depth","variant_quality","MitoTip_score","MitoTip_interpretation","QUAL",
"FILTER","MQM_INFO","MQMR_INFO", "QA_INFO", "QR_INFO", "SAF_INFO", "SAR_INFO","SRF_INFO","SRR_INFO", 
"SBR_INFO","SBA_INFO","POS_FILTER","SBR_FILTER","SBA_FILTER", "MQMR_FILTER","AQR_FILTER","GT_FORMAT",
"QR_FORMAT","AQR_FORMAT","QA_FORMAT","AQA_FORMAT","INFO","FORMAT"]

#remove columns and store the remaining cols in a new df called final_report
final_report=file_df.drop(remove_cols,axis=1) 

#Concatenate the CV_gnomad_annot dataframe to final_report dataframe
final_report=pd.concat([final_report,CV_gnomad_annot],axis=1)

#Re-arrange columns of the final_report
final_report=final_report[['CHR','POS','REF','ALT','SAMPLE', 'HGVS', 'gene/locus', 'gene/locus description', 'variant_heteroplasmy', 
'alt_depth', 'ref_depth', 'total_locus_depth', 'COHORT_COUNT', 'tier', 'commercial_panels', 'phylotree_haplotype', 
'gnomad_mitotip_score', 'gnomad_mitotip_trna_prediction', 'MitoTip_percentile', 'anticodon', 'allele_frequency_mitomap', 
'disease_mitomap', 'MGRB_frequency', 'MGRB_FILTER', 'MGRB_AC', 'MGRB_AN', 'phylotree_mut', 'locus_mitomap', 
'num_references_mitomap', 'variant_amino_acid_change_mitomap', 'codon_position_mitomap', 'codon_number_mitomap', 
'num_disease_references_mitomap', 'RNA_mitomap', 'homoplasmy_mitomap', 'heteroplasmy_mitomap', 'status_mitomap', 
'disease_amino_acid_change_mitomap', 'GenBank_frequency_mitomap','clinvar_pathogenic', 'gnomad_mitotip_score', 
'gnomad_mitotip_trna_prediction', 'AC_hom', 'AC_het', 'AF_hom', 'AF_het', 'max_hl']]

print(final_report)

#write the final report as a csv file
final_report.to_csv("1760_mity_vcfanno_final_report.csv",index=False)