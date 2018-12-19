.separator “\t”

create table estr_gtex(
chrom text not null,
start integer not null,
end integer not null,
str_id text not null,
gene_id text,
gene_name text,
tissue text,
beta real,
pval real,
caviar real);

.import /storage/ryanicky/webstr_files/estr_gtex_webstr.tab estr_gtex



create table mutrates(
chrom text not null,
start integer not null,
end integer not null,
str_id text not null,
est_logmu_ml real,
est_beta_eff_ml real,
est_beta_ml real,
est_pgeom_ml real,
stderr_ml real,
numsamples_ml integer,
up real,
down real,
p real,
period integer,
motif text,
uninterrupted_length real,
pred_mu_1 real,
pred_mu_se_1 real,
zscore_1 real,
pred_mu_2 real,
pred_mu_se_2 real,
zscore_2 real,
filter1 text);


.import /storage/ryanicky/webstr_files/mutrates_webstr.tab mutrates



create table allelstattemp(
locus1 text,
locus2 text,
allele integer,
freq_het real,
MAF real,
KL real,
r2 real,
pval real);

.import /storage/ryanicky/webstr_files/str_imputation_allelestats.tab allelstattemp


create table locstattemp(
ctr integer,
pos text not null,
str_id text not null,
loo_r real,
loo_pVal real,
loo_concordance real,
loo_numSamples real,
wgs_eur_pVal real,
wgs_eur_r real,
wgs_eur_concordance real,
wgs_eur_numSamples real,
wgs_afr_pVal real,
wgs_afr_r real,
wgs_afr_concordance real,
wgs_afr_numSamples real,
wgs_eas_pVal real,
wgs_eas_r real,
wgs_eas_concordance real,
wgs_eas_numSamples real,
affy_eur_pVal real,
affy_eur_r real,
affy_eur_concordance real,
affy_eur_numSamples real,
affy_afr_pVal real,
affy_afr_r real,
affy_afr_concordance real,
affy_afr_numSamples real,
affy_eas_pVal real,
affy_eas_r real,
affy_eas_concordance real,
affy_eas_numSamples real,
omni_eur_pVal real,
omni_eur_r real,
omni_eur_concordance real,
omni_eur_numSamples real,
omni_eas_pVal real,
omni_eas_r real,
omni_eas_concordance real,
omni_eas_numSamples real);

.separator “,”

.import /storage/ryanicky/webstr_files/str_imputation_locstats.csv locstattemp



CREATE TABLE allelstat AS
       select lt.str_id, a.*
       FROM locstattemp lt, allelstattemp a
       WHERE A.locus1 = lt.pos;


CREATE TABLE locstat AS
select pos,str_id,loo_r,loo_concordance,wgs_eur_r,wgs_eur_concordance,wgs_afr_r,wgs_afr_concordance,wgs_eas_r,wgs_eas_concordance from locstattemp;


drop table allelstattemp;
drop table locstattemp;

create index allelstat_strid_idx on allelstat(str_id);
create index locstat_strid_idx on locstat(str_id);
create index mutrates_strid_idx on mutrates(str_id);
create index estr_gtex_strid_idx on estr_gtex(str_id);

