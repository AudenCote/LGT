args = commandArgs(trailingOnly = TRUE)

lib_path <- "/Library/Frameworks/R.framework/Versions/3.5/Resources/library"

library(rvest, lib.loc = lib_path)
library(httr, lib.loc = lib_path)

baseurl <- "http://revigo.irb.hr/"

goList = args[1]

revigo_session <- html_session("http://revigo.irb.hr/")
revigo_form <- html_form(revigo_session)[[1]]  
filled_form <- set_values(revigo_form,'goList' = goList, 'cutoff'=0.4, 'isPValue'="no")
result_page <- submit_form(revigo_session, filled_form, submit='startRevigo')

results_table <- html_table(result_page)[[1]]
names(results_table) <- results_table[2,]
results_table <- results_table[3:nrow(results_table),]

write.csv(results_table, 'revigo_results')
