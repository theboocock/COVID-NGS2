gisaid_submission = c("submitter"	,
  "fn",
  "covv_virus_name",
  "covv_type",
  "covv_passage",
  "covv_collection_date",
  "covv_location",
  "covv_add_location",
  "covv_host",
  "covv_add_host_info",
  "covv_gender",
  "covv_patient_age",
  "covv_patient_status",
  "covv_specimen",
  "covv_outbreak",
  "covv_last_vaccinated",
  "covv_treatment",
  "covv_seq_technology",
  "covv_assembly_method",
  "covv_coverage",
  "covv_orig_lab",
  "covv_orig_lab_addr",
  "covv_provider_sample_id",
  "covv_subm_lab",
  "covv_subm_lab_addr",
  "covv_subm_sample_id",
  "covv_authors")

username = "theboocock"
fasta = "all_sequences.fasta"
covv_type = "betacoronavirus"
location_string = "North America / USA / California / Los Angeles County"
gisaid_submit = gisaid_phased_filt %>% filter(submitting_lab == "UCLA") %>% mutate(strain_new = str_replace_all(strain,"LA-[0]+",""))  %>% mutate(name_out=paste("hCoV-19/USA/CA-UCLA",strain_new,"/2020",sep=""))
fasta_in_la = readDNAStringSet("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/outputs/final/merged/all.fasta")
fasta_in_la = fasta_in_la[gisaid_submit$strain]

#### read 

lineages_remake_fasta = read.csv("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/lineages.csv")

names(fasta_in_la) = gisaid_submit$name_out
#gisaid_submit$strain

fasta_in_la2 = readDNAStringSet("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/fastas/all.fasta")
fasta_in_la2 =fasta_in_la2[gisaid_submit$strain]
names(fasta_in_la2) =  gisaid_submit$name_out
gisaid_sub_df = gisaid_submit %>% mutate(submitter=username,fn=fasta,covv_virus_name=name_out, 
                         covv_passage="Original",covv_collection_date=date_fix,
                         covv_location="North America / USA / California / Los Angeles County",
                         covv_host = "Human",covv_gender="unknown",
                         covv_patient_age="unknown",covv_patient_status="unknown",
                         covv_seq_technology="Illumina",
                         covv_orig_lab="UCLA Pathology Clinical Microbiology Lab",
                         covv_orig_lab_addr="Brentwood Annex,11633 San Vicente Blvd, Los Angeles, CA 90049",
                         covv_subm_lab="Kruglyak Lab",covv_subm_lab_addr="695 Charles E Young Drive South, Los Angeles, CA 90095, US",covv_authors="Guo et al.",covv_treatment="",
                         covv_last_vaccinated="",covv_subm_sample_id=strain,covv_add_host_info="",covv_assembly_method="",covv_add_location="",
                         covv_provider_sample_id="",covv_outbreak="",covv_coverage="",covv_specimen="",covv_type=covv_type) %>% select(one_of(gisaid_submission))





writeXStringSet(fasta_in_la,filepath = "data/all_sequences.fasta")
writeXStringSet(fasta_in_la2,filepath = "data/all_sequences2.fasta")

write.csv(gisaid_sub_df,file="data/gisaid_sub_df.csv", quote=T,row.names=F)

g = gisaid_phased_filt %>% filter(submitting_lab != "UCLA" & location == "Los Angeles County")
