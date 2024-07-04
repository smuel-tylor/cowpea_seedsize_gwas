library(tidyverse)
library(patchwork)
library(ComplexUpset)
#library(UpSetR)

read_blink <- \(fp){
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
probit_meta <- \(gwas_list){
  
  z_squared <- \(p_value){
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

gwas_blink_list_mta <- gwas_blink_list |>
  map(\(x) arrange(x, Chromosome, Position) |>
        mutate(Position_G = cumsum(Position),
               #this is slightly larger than 5.92
               mta = ifelse(neg_log10p > -log(0.05 / n_bfc, base = 10),
                            "sig", "ns"),
        )
  )

gwas_blink_list_mta |>
  map(\(x) x|>
        ggplot(aes(Position_G, neg_log10p, colour = Chromosome)) + 
        geom_point() +
        ylim(0, 25)
  )

gwas_blink_list_sig <- gwas_blink_list_mta |>
  map(\(x) filter(x, mta == "sig"))
#should mean we now have different length sets in each element of the list
map_vec(gwas_blink_list_sig, nrow)

#function that scans along sig column to id changes
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

gwas_blink_list_psnp <- gwas_blink_list_mta |>
  map(last_threshold) |>
  map(\(x) x |>
        filter(mta == "sig") |>
        group_by(last_threshold) |>
        mutate(max_neg_log10p = max(neg_log10p),
               psnp = ifelse(max_neg_log10p == neg_log10p,
                             "psnp",
                             "snp"
               )
        ) |>
        ungroup() |>
        filter(psnp == "psnp") |>
        #this last transformation is needed to ensure upset works
        as.data.frame()
  )

gwas_blink_list_sig_snps <- map(gwas_blink_list_sig, \(x) x[ , "SNP"])

ur <- UpSetR::fromList(gwas_blink_list_sig_snps)
upset(ur,
      intersect = names(ur)
      )
#not working - is it an R update issue?

gwas_blink_list_psnp_snps <- map(gwas_blink_list_psnp, \(x) x[ , "SNP"])
upset(fromList(gwas_blink_list_psnp_snps),
      list(names(fromList(gwas_blink_list_psnp_snps)))
)

#MLM
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
               #this is slightly larger than 5.92
               mta = ifelse(neg_log10p > -log(0.05 / n_bfc, base = 10),
                            "sig", "ns"),
               )
  )

gwas_mlm_list_mta |>
        map(\(x) x|>
            ggplot(aes(Position_G, neg_log10p, colour = Chromosome)) + 
        geom_point() +
        ylim(0, 25)
        )

#This confirms that old psnp analysis was inappropriate
#11 p-snps should be present in just this segment of chr3
#The previous analysis showed only 13 across the entire genome
gwas_mlm_list_mta$MLM_meta |>
  filter(Chromosome == "VU03") |>
  ggplot(aes(Position_G, neg_log10p, colour = Chromosome)) + 
  geom_line() +
  ylim(0, 25) +
  xlim(1.4775e11, 1.483e11)

gwas_mlm_list_sig <- gwas_mlm_list_mta |>
  map(\(x) filter(x, mta == "sig"))
#should mean we now have different length sets in each element of the list
map_vec(gwas_mlm_list_sig, nrow)

#this segment won't handle the missing parts of the list
gwas_mlm_list_psnp <- gwas_mlm_list_mta[map_vec(gwas_mlm_list_sig, nrow) != 0] |>
  map(last_threshold) |>
  map(\(x) x |>
        filter(mta == "sig") |>
        group_by(last_threshold) |>
        mutate(max_neg_log10p = max(neg_log10p, na.rm = TRUE),
               psnp = ifelse(max_neg_log10p == neg_log10p,
                             "psnp",
                             "snp"
               )
        ) |>
        ungroup() |>
        filter(psnp == "psnp") |>
        as.data.frame()
  )

#should mean we now have different length sets in each element of the list
map_vec(gwas_mlm_list_psnp, nrow)

gwas_mlm_list_sig_snps <- map(gwas_mlm_list_sig, \(x) x[ , "SNP"])
upset(fromList(gwas_mlm_list_sig_snps),
      list(names(fromList(gwas_mlm_list_sig_snps)))
)

gwas_mlm_list_psnp_snps <- map(gwas_mlm_list_psnp, \(x) x[ , "SNP"])
#padding this
gwas_mlm_list_psnp_snps$MLM_field <- character(0)
gwas_mlm_list_psnp_snps$MLM_width <- character(0)
gwas_mlm_list_psnp_snps$MLM_length <- character(0)
upset(fromList(gwas_mlm_list_psnp_snps),
      list(names(fromList(gwas_mlm_list_psnp_snps)))
)

mlm_bnk_sig <- c(gwas_mlm_list_sig_snps, gwas_blink_list_sig_snps)
names(mlm_bnk_sig) <- str_replace_all(names(mlm_bnk_sig), "_", " ")
a <- upset(fromList(mlm_bnk_sig),
           intersect = c("BLINK density", "BLINK length", "BLINK width",
                         "BLINK cvars", "BLINK field", "BLINK gh", "BLINK meta",
                         "MLM density", "MLM length", "MLM width",
                         "MLM cvars", "MLM field", "MLM gh", "MLM meta"
           ),
           sort_sets = FALSE
           #need to re-write the below consistent with ComplexUpset
      #nsets = length(mlm_bnk_sig),
      #sets = names(mlm_bnk_sig)[length(mlm_bnk_sig):1],
      #keep.order = TRUE,
      #mainbar.y.max = 135,
      #order.by = c("freq"),
      #sets.bar.color = rep(c(gray(0.2), gray(0.8)), each = 7),
      #set_size.scale_max = 170
)

a

#need to do psnps
mlm_bnk_psnp <- c(gwas_mlm_list_psnp_snps, gwas_blink_list_psnp_snps)

mlm_bnk_psnp_table <- list_rbind(
  list(
    list_rbind(gwas_blink_list_psnp, names_to = "analysis"),
    list_rbind(gwas_mlm_list_psnp, names_to = "analysis")
  )
)
write.csv(mlm_bnk_psnp_table, "mlm_blink_psnp.csv")

distinct(mlm_bnk_psnp_table, SNP, Chromosome, Position) |>
  arrange(Chromosome, Position)

names(mlm_bnk_psnp) <- str_replace_all(names(mlm_bnk_psnp), "_", " ")
b <- upset(upsetr::fromList(mlm_bnk_psnp),
           nsets = length(mlm_bnk_psnp),
           #sets = names(mlm_bnk_psnp)[length(mlm_bnk_psnp):1],
           sets = c("MLM density", "MLM length", "MLM width",
                    "MLM cvars", "MLM field",
                    "MLM gh", "MLM meta",
                    "BLINK density", "BLINK length", "BLINK width",
                    "BLINK cvars", "BLINK field", "BLINK gh", "BLINK meta"
           ),
           keep.order = TRUE,
           mainbar.y.max = 135,
           order.by = c("freq"), decreasing = FALSE,
           sets.bar.color = rep(c(gray(0.2), gray(0.8)), each = 7),
           set_size.show = TRUE,
           set_size.scale_max = 165,
           set_sizes = (upset_set_size() + xlim(-20, 165))
)

b

pdf("blink_mlm_upset.pdf", paper = "a4")

(
(plot_spacer() + wrap_elements(a$Main_bar) + plot_layout(widths = c(1, 3))) /
    (
      wrap_elements(a$Sizes) +
       (wrap_elements(a$Matrix) / plot_spacer() + plot_layout(heights = c(10, 1))) +
       plot_layout(widths = c(1, 3))
     )
) /
  (
    (plot_spacer() + wrap_elements(b$Main_bar) + plot_layout(widths = c(1, 3))) /
  (
    wrap_elements(b$Sizes) +
     (wrap_elements(b$Matrix) / plot_spacer() + plot_layout(heights = c(10, 1))) +
    plot_layout(widths = c(6, 3))
  )
  ) +
  plot_layout(heights = c(1, 1, 2))

dev.off()  
