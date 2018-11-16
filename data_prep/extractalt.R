library(dplyr)
library(iotools)
library(tidyr)


#### Arg 1 should be filename of input text with full path.
#### arg 2 should be output filename with full path.

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} 

print(args[1])

dat1 <- read.csv.raw(file=as.character(args[1]),sep="|",header=F)
sum1 <- dat1 %>% select(V1,V2,V3,V4,V5,V6) %>% group_by(V2,V3,V4,V5,V6) %>% summarize(cnt=n())

t1 <- sum1 %>% mutate(V6 =strsplit(V6, ",")) %>% unnest(V6)
t1df <- data.frame(t1)
t1df$V6 <- trimws(t1df$V6)
t2df = t1df %>% group_by (V2, idx = cumsum(V4 == 1L)) %>% mutate(alt_num = row_number()) %>% ungroup %>% select (-idx)

names(t2df) <- c("POS","ID","REF","FILTER","cnt","ALT","ALT_NUM")
write.table(t2df[,c("POS","ID","REF","FILTER","ALT","ALT_NUM")],sep=",",row.names=F,file=args[2],quote=F,col.names=F)
