

PRAGMA TEMP_STORE=2;

CREATE TABLE altGT AS
select * from
(select sample,str_id,altref_gt1 as altref_gt
from tempGT
UNION
select sample,str_id,altref_gt2 as altref_gt
from tempGT
);


create index altGT_strid_idx on altGT(str_id);

create index homozyg_strid_idx on vcfhomozyg(str_id);

drop table tempGT; 

create index vcfAlt_strid_idx on vcfAlt(str_id);

create table vcfBase(
pos integer,
str_id text,
ref text);

INSERT INTO vcfBase SELECT pos,str_id,ref from altemp group by pos,str_id,ref;

create index vcfBase_strid_idx on vcfBase(str_id);


vacuum;

