library(tidyverse);
library(RColorBrewer);
if(!"make_plots"%in%ls()) {make_plots = F;}

#Case 23 ----
#Read in data; remove tail columns that are empty
case23_ranked <- read_csv(file="case23_20Dec2014.csv",col_names=F,na="",skip=4,n_max=32);
colnames(case23_ranked) <- colnames(read_csv(file="case23_20Dec2014.csv",col_names=T,na="",skip=2,n_max=1));
case23_ranked <- select(case23_ranked, -starts_with("X"));
colnames(case23_ranked) <- 
  colnames(case23_ranked) %>%
  str_to_upper() %>%
  str_replace_all(pattern = fixed("."),replacement = " ");

#Recharacterize as ordered
case23_ordered = t(apply(case23_ranked,1,order));
for(i in 1:nrow(case23_ordered)) {
  case23_ordered[i,(sum(!is.na(case23_ranked[i,]))+1):ncol(case23_ranked)] = NA;
}
case23_ordered = case23_ordered[,which(colSums(!apply(case23_ordered,2,is.na))>0)];
if(make_plots) {
  ordered_barplot(case23_ranked, 
                  col_palette_range = c("#800026", "#FFCB05"), 
                  file_path = "../../barplot_case23.png",
                  title_text = "Case A");
}

#Case 111 ----
#Read in data; remove tail columns that are empty
case111_ranked <- read_csv(file="case111_20Dec2014.csv",col_names=F,na="",skip=4,n_max=32);
colnames(case111_ranked) <- colnames(read_csv(file="case111_20Dec2014.csv",col_names=T,na="",skip=2,n_max=1));
case111_ranked <- select(case111_ranked, -starts_with("X"));
colnames(case111_ranked) <- 
  colnames(case111_ranked) %>%
  str_to_upper() %>%
  str_replace_all(pattern = fixed("."),replacement = " ");

#Recharacterize as ordered
case111_ordered = t(apply(case111_ranked,1,order));
for(i in 1:nrow(case111_ordered)) {
  case111_ordered[i,(sum(!is.na(case111_ranked[i,]))+1):ncol(case111_ranked)] = NA;
}
case111_ordered = case111_ordered[,which(colSums(!apply(case111_ordered,2,is.na))>0)];
if(make_plots) {
  ordered_barplot(case111_ranked, 
                  col_palette_range = c("#00441B", "#FFCB05"), 
                  file_path = "../../barplot_case111.png",
                  title_text = "Case B");
}

#Case 83 ----
#Read in data; remove tail columns that are empty
case83_ranked <- read_csv(file="case83_20Dec2014.csv",col_names=F,na="",skip=4,n_max=32);
colnames(case83_ranked) <- colnames(read_csv(file="case83_20Dec2014.csv",col_names=T,na="",skip=2,n_max=1));
case83_ranked <- select(case83_ranked, -starts_with("X"));
colnames(case83_ranked) <- 
  colnames(case83_ranked) %>%
  str_to_upper() %>%
  str_replace_all(pattern = fixed("."),replacement = " ");

#Recharacterize as ordered
case83_ordered = t(apply(case83_ranked,1,order));
for(i in 1:nrow(case83_ordered)) {
  case83_ordered[i,(sum(!is.na(case83_ranked[i,]))+1):ncol(case83_ranked)] = NA;
}
case83_ordered = case83_ordered[,which(colSums(!apply(case83_ordered,2,is.na))>0)];
if(make_plots) {
  ordered_barplot(case83_ranked, 
                  col_palette_range = c("#00274C", "#FFCB05"), 
                  file_path = "../../barplot_case83.png",
                  title_text = "Case C");
}

