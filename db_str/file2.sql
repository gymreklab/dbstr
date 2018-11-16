

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
ref text,
filter text);

INSERT INTO vcfBase SELECT pos,str_id,ref,filter from vcfAlt group by pos,str_id,ref,filter;

create index vcfBase_strid_idx on vcfBase(str_id);

CREATE VIEW Sample_szes  AS
select base.*,hzyg.sample,' ' alt,length(base.ref) lent,2 mult from
            vcfBase as base,
            vcfhomozyg as hzyg
            where
                base.str_id=hzyg.str_id

            UNION ALL

            select base.*,GT.sample,alt.alt,length(alt.alt) lent,1 mult from
            vcfBase as base,
            vcfAlt as alt,
            altGT as gt
            where
                base.pos=alt.pos
            and base.str_id=alt.str_id
            and base.str_id = gt.str_id
            and alt.str_id = gt.str_id
            and gt.altref_gt = alt.altorder;



vacuum;

