
rm(list = ls())
require("dplyr")
require("tidyr")
require("reshape2")
require("data.table")
require("optparse")

option_list <- list(make_option(c("-i", "--input"), type="character", default=NULL, help="filtered_cnv_input", metavar="character"),
                    make_option(c("-o", "--outfile"), type="character", default=NULL, help="tip_probability_output", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#print(names(opt))
#print(opt)

# opt$input = "../data/cadd456b3a9334a5/medicc_input.csv"
# opt$outfile = "../data/cadd456b3a9334a5/lazac_input.csv"

read_input = function(input) {
  df = read.delim(input) %>% 
    dplyr::mutate(bin = paste(chrom, start, end, sep = "_"))
  
  df %>% 
    dplyr::select(sample_id, cn_a, cn_b, chrom, start, end) %>% 
    dplyr::rename(node = sample_id)
}

print(paste0("Reading ", opt$input))
print(paste0("Saving to ", opt$outfile))
df = read_input(input = opt$input)
data.table::fwrite(df, opt$outfile, row.names = F, quote = F)
