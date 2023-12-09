library(UpSetR)
library(tidyverse)

read_blink <- function(fp){
  read.csv(fp) |>
    select(SNP:P.value) |>
    mutate(neg_log10p = -log(P.value, base = 10),
           Trait = str_to_lower(
             str_extract(fp, "Blink.([A-Za-z]+).GWAS", group = 1)
           ),
           Trait = case_match(Trait,
             "Seeddensity" ~ "density",
             "Average10seedslength" ~ "length",
             "SeedWeightCVARS" ~ "cvars",
             "SeedWeightField" ~ "field",
             "SeedWeightGH" ~ "gh",
             "Average10seedswidth" ~ "width"
           )
           )
}

blink_dirs <- str_c("copy_of_MWK_dir/GWAS_Results/Seed ",
                  c("Density", "Length",
                    "weight CVARS", "weight Field", "weight GH",
                    "width"
                    ),
                  "/GAPIT"
                  )

blink_fp <- map_chr(blink_dirs, \(x) list.files(x, full.names = TRUE))

gwas_blink_list <- map(blink_fp, read_blink)
names(gwas_blink_list) <- str_c(
  "BLINK ",
  c("density", "length", "cvars", "field", "gh", "width")
)
  
#what is n for bonferroni
n_bfc <- unique(map_vec(gwas_blink_list, nrow))

gwas_blink_list[[1]] |>
  arrange(Chromosome, Position) |>
  mutate(Position_G = cumsum(Position)) |>
ggplot(aes(Position_G, neg_log10p, colour = Chromosome)) + 
  geom_point()

gwas_blink_list_sig <- gwas_blink_list |>
  map(\(x) filter(x, neg_log10p > -log(0.05 / n_bfc, base = 10)))
#should mean we now have different length sets in each element of the list
map_vec(gwas_blink_list_sig, nrow)

gwas_blink_list_snps <- map(gwas_blink_list_sig, \(x) x[ , "SNP"])
upset(fromList(gwas_blink_list_snps))

read_tassel <- function(fp){
  read.delim(fp) |>
    filter(!is.na(Pos)) |> #removes a null row
    select(Trait:Pos, p) |>
    rename(SNP = Marker,
           Chromosome = Chr,
           Position = Pos,
           P.value = p
           ) |>
    mutate(neg_log10p = -log(P.value, base = 10),
           Position = as.numeric(Position)
           )
}

mlm_dirs <- str_c("copy_of_MWK_dir/GWAS_Results/Seed ",
                    c("Density", "Length",
                      "weight CVARS", "weight Field", "weight GH",
                      "width"
                    ),
                    "/Tassel"
)

mlm_fp <- map_chr(mlm_dirs, list.files, full.names = TRUE)

gwas_mlm_list <- map(mlm_fp, read_tassel)
#what is n for bonferroni
unique(map_vec(gwas_mlm_list, nrow)) == n_bfc

gwas_mlm_list[[1]] |>
  arrange(Chromosome, Position) |>
  mutate(Position_G = cumsum(Position)) |>
  ggplot(aes(Position_G, neg_log10p, colour = Chromosome)) + 
  geom_point()

gwas_mlm_list_sig <- gwas_mlm_list |>
  map(\(x) filter(x, neg_log10p > -log(0.05 / n_bfc, base = 10)))
#should mean we now have different length sets in each element of the list
map_vec(gwas_mlm_list_sig, nrow)

a <- map(gwas_blink_list_sig, \(x) x[ , "SNP"])
names(a) <- c("A", "B", "C", "D", "E", "F")
upset(fromList(a))
