library(dplyr)
library(data.table)
options(stringsAsFactors = F)

dir = '/scratch/t.cri.yliang/ptrs-ukb/prs_on_cri'
output = '/gpfs/data/im-lab/nas40t2/yanyul/PTRS/prs_cri_subset1.tsv.gz'

traits = read.table('trait_list.txt', header = F)$V1
ranges = read.table('range_file.txt', header = F)$V1

indivs = NULL
dd_all = NULL
for(trait in traits) {
  for(range in ranges) {
    ddr = NULL
    for(chr in 1 : 22) {
      message(glue::glue('Working on {trait} {range} chr{chr}'))
      fn = glue::glue('{dir}/{trait}.chr{chr}.{range}.sscore')
      tmp = fread(fn, header = T)
      coln = glue::glue('{trait}_x_{range}')
      colnames(tmp) = c('indiv', coln)
      if(is.null(ddr)) {
        ddr = tmp
        indivs = tmp$indiv
      } else {
        tmp = left_join(data.frame(indiv = indivs), tmp, by = 'indiv')
        if(sum(is.na(tmp[[coln]])) > 0) {
          message('Not the same set of individuals.')
          quit()
        }
        ddr[[coln]] = ddr[coln] + tmp[[coln]]
      }
    }
    if(is.null(dd_all)) {
      dd_all = ddr
    } else {
      dd_all = inner_join(dd_all, ddr, by = 'indiv')
    }
  }
}
dd_all = dd_all %>% rename(s = indiv)
gz <- gzfile(output, 'w')
write.table(dd_all, gz, quote = F, col = T, row = F, sep = '\t')
close(gz)
