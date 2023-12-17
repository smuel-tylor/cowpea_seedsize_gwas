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
  "BLINK_",
  c("density", "length", "cvars", "field", "gh", "width")
)

#taking as input the list components relevant to the meta-analysis
#this function carries out the probit transformation
probit_meta <- function(gwas_list){
  
  z_squared <- function(p_value){
    qnorm(1 - 0.5 * p_value)^2
  }
  
  imap(gwas_list,
       \(x, idx) rename_with(x,
                             ~ str_c(idx, "_", .x),
                             !contains(c("SNP", "Chromosome", "Position"))
       )
  ) |>
    reduce(full_join) |>
    select(contains(c("SNP", "Chromosome", "Position", "P.value", "log10p"))) |>
    arrange(Chromosome, Position) |>
    mutate(across(contains("P.value"), z_squared, .names = "{.col}_z_sq"),
           meta_Chi_sq = rowSums(across(contains("z_sq"))),
           meta_P.value = pchisq(meta_Chi_sq,
                                        length(gwas_list),
                                        lower.tail = FALSE
           ),
           meta_neg_log10p = -log(meta_P.value, base = 10)
    )
}

gwas_blink_list$BLINK_meta <- probit_meta(
  gwas_blink_list[c("BLINK_cvars", "BLINK_field", "BLINK_gh")]
  ) |>
  select(c("SNP", "Chromosome", "Position", "meta_P.value", "meta_neg_log10p")
         ) |>
  rename(P.value = meta_P.value,
         neg_log10p = meta_neg_log10p
         )

#what is n for bonferroni
n_bfc <- unique(map_vec(gwas_blink_list, nrow))

gwas_blink_list |>
map(\(x) arrange(x, Chromosome, Position) |>
  mutate(Position_G = cumsum(Position)) |>
ggplot(aes(Position_G, neg_log10p, colour = Chromosome)) + 
  geom_point() +
  ylim(0, 25)
)
#There's an Inf in the mix from the meta-analysis...
gwas_blink_list$BLINK_meta |> filter(neg_log10p > 15)

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

gwas_mlm_list <- map(mlm_fp, \(x) read_tassel(x) |>
                       select(c("SNP", "Chromosome", "Position",
                                "P.value", "neg_log10p"
                                )
                              )
)
names(gwas_mlm_list) <- str_c(
  "MLM_",
  c("density", "length", "cvars", "field", "gh", "width")
)

gwas_mlm_list$MLM_meta <- probit_meta(
  gwas_mlm_list[c("MLM_cvars", "MLM_field", "MLM_gh")]
) |>
  select(c("SNP", "Chromosome", "Position", "meta_P.value", "meta_neg_log10p")
  ) |>
  rename(P.value = meta_P.value,
         neg_log10p = meta_neg_log10p
  )

#what is n for bonferroni
unique(map_vec(gwas_mlm_list, nrow)) == n_bfc

gwas_mlm_list_mta <- gwas_mlm_list |>
  map(\(x) arrange(x, Chromosome, Position) |>
        mutate(Position_G = cumsum(Position),
               mta = ifelse(neg_log10p > -log(0.05 / n_bfc, base = 10),
                            "sig", "ns"),
               )
  )

gwas_mlm_list_mta |>
        map(\(x) ggplot(x, aes(Position_G, neg_log10p, colour = Chromosome)) + 
        geom_point() +
        ylim(0, 25)
        )

gwas_mlm_list_sig <- gwas_mlm_list_mta |>
  map(\(x) filter(x, mta == "sig"))
#should mean we now have different length sets in each element of the list
map_vec(gwas_mlm_list_sig, nrow)

#function that scans along sig column to id changes
####wip - struggling with it a bit
last_threshold <- \(x){
  x[ , "last_threshold"] <- rep(NA, nrow(x))
  for (i in 1:nrow(x)){
    x[i, "last_threshold"] <- ifelse(i == 1,
                                     x[i, "SNP"],
                                     ifelse(
                                       x[i - 1, "mta"] == x[i, "mta"],
                                       x[i - 1, "last_threshold"],
                                       x[i, "SNP"]
                                     )
    )
  }
  x
}

#couldn't work out how to do this using modify... kept throwing errors
  # #last <-
  #   modify2(rep(NA, nrow(x)), 1:length(last),
  #         \(y, i) y[i] <- ifelse(i == 1,
  #                                x[i, "SNP"],
  #                                ifelse(
  #                                    x[i - 1, "mta"] == x[i, "mta"],
  #                                    y[i - 1],
  #                                    x[i, "SNP"]
  #                                  )
  # )
  # )
  #data.frame(x, last_threshold = last)

#below continues to be an issue
gwas_mlm_list_psnp <- gwas_mlm_list_mta |>
  map(last_threshold) |>
  map(\(x) split(x, ~ last_threshold) |>
        map(\(y) mutate(y, psnp = ifelse(max(neg_log10p), "psnp", "snp"))
            ) |>
        unsplit(~ last_threshold)
  ) |>
  map(\(x) filter(x, mta == "sig"))


        filter(max_neg_log10p == "psnp")
#should mean we now have different length sets in each element of the list
map_vec(gwas_mlm_list_mta, nrow)

gwas_mlm_list_snps <- map(gwas_mlm_list_mta, \(x) x[ , "SNP"])
upset(fromList(gwas_mlm_list_snps))

mlm_bnk <- c(gwas_mlm_list_snps, gwas_blink_list_snps)
names(mlm_bnk) <- str_replace_all(names(mlm_bnk), "_", " ")
a <- upset(fromList(mlm_bnk),
      nsets = length(mlm_bnk),
      #sets = names(mlm_bnk)[length(mlm_bnk):1],
      sets = c("BLINK density", "BLINK length", "BLINK width",
               "BLINK cvars", "BLINK field", "BLINK gh", "BLINK meta",
               "MLM density", "MLM length", "MLM width",
               "MLM cvars", "MLM field", "MLM gh", "MLM meta"
               ),
      keep.order = TRUE,
      order.by = c("freq"),
      sets.bar.color = rep(c(gray(0.2), gray(0.8)), each = 7)
)

a

#need to do psnps
gwas_mlm_list_sig
#need to do something involving lags to establish distances...
  