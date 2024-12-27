library(ComplexUpset)
library(tidyverse)
library(patchwork)

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
    ) |>
    arrange(Chromosome, Position)
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
           #added a fix here because this was calculating p = 0
           meta_neg_log10p = ifelse(meta_P.value != 0,
                                    -log(meta_P.value, base = 10),
                                    -log(2.22e-16, base = 10)
           )
    )
}

#And now we tack on the meta-analysis result to the rest of the GWAS outcomes
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
#generate absolute positions on genome
Chr_max <- gwas_blink_list[[1]] |>
  arrange(Chromosome, Position) |>
  group_by(Chromosome) |>
  mutate(Chromosome_max = max(Position)) |>
  ungroup() |>
  distinct(Chromosome, Chromosome_max) |>
  mutate(Chromosome_max_Genome = cumsum(Chromosome_max))

# gwas_blink_list |>
# map(\(x) full_join(x, Chr_max) |>
#       arrange(Chromosome, Position) |>
#       mutate(Position_G = Chromosome_max_Genome - (Chromosome_max - Position)
#              ) |>

# gwas_blink_list |>
#   map(\(x) arrange(x, Chromosome, Position) |>
#         mutate(Position_G = cumsum(Position)) |>
#         ggplot(aes(Position_G, neg_log10p, colour = Chromosome)) + 
#         geom_point() +
#         ylim(0, 25)
#   )

gwas_blink_list_mta <- gwas_blink_list |>
  map(\(x) full_join(x, Chr_max) |>
        arrange(Chromosome, Position) |>
        mutate(Position_G = Chromosome_max_Genome - (Chromosome_max - Position),
               #this is slightly larger than 5.92
               mta = ifelse(neg_log10p > -log(0.05 / n_bfc, base = 10),
                            "sig", "ns"),
        )
  )

#write a file with all the SNPs listed
SNPs <- gwas_blink_list_mta[[1]] |>
  select(SNP:Position, Position_G) |>
  arrange(Position_G)

write.csv(SNPs, "SNPs.csv", row.names = FALSE)

#plot out Manhattans
gwas_blink_list_mta |>
  map(\(x) x|>
        ggplot(aes(Position_G, neg_log10p, colour = Chromosome)) +
        geom_point() +
        ylim(0, 16)
      )

gwas_blink_list_sig <- gwas_blink_list_mta |>
  map(\(x) filter(x, mta == "sig"))
#should mean we now have different length sets in each element of the list
map_vec(gwas_blink_list_sig, nrow)

#function that scans along sig column to id changes
#1224, added failsafe arrange to make sure ordering is correct
last_threshold <- \(x){
  x <- arrange(x, Position_G)
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
  return(x)
}

#function that then calls psnps based on this scan
call_psnp <- \(x) x |>
  filter(mta == "sig") |>
  group_by(last_threshold) |>
  #Here I'm making sure this value doesn't set to -Inf
  # so that rbind will work later
  mutate(max_neg_log10p = ifelse(max(neg_log10p) == -Inf, -1, max(neg_log10p)),
         psnp = ifelse(max_neg_log10p == neg_log10p,
                       "psnp",
                       "snp"
         )
  ) |>
  ungroup() |>
  filter(psnp == "psnp") |>
  mutate(psnp = as.character(psnp))

#use these functions and generate a tibble
gwas_blink_list_psnp <- gwas_blink_list_mta |>
  map(\(x) split(x, factor(x$Chromosome))) |>
  map_depth(2, last_threshold) |>
  map_depth(2, call_psnp) |>
  map_depth(1, list_rbind)

#reformat and use complexupset to plot
gwas_blink_sig_snps <- gwas_blink_list_sig |>
  list_rbind(names_to = "method_phenotype") |>
  select(-c(P.value:Trait)) |>
  pivot_wider(names_from = method_phenotype, values_from = mta) |>
  mutate(across(starts_with("BLINK"),
                ~ ifelse(.x == "sig", 1, 0) |> replace_na(0)
  )
  )

upset(gwas_blink_sig_snps,
      intersect = c("BLINK_density", "BLINK_length", "BLINK_width",
                    "BLINK_cvars", "BLINK_field", "BLINK_gh", "BLINK_meta"
      )
)

gwas_blink_psnp_snps <- gwas_blink_list_psnp |>
  map_depth(1, \(x) arrange(x, Position_G)) |>
  list_rbind(names_to = "method_phenotype") |>
  select(-c(P.value:Trait, mta:max_neg_log10p)) |>
  pivot_wider(names_from = method_phenotype, values_from = psnp) |>
  mutate(across(starts_with("BLINK"),
                ~ ifelse(.x == "psnp", 1, 0) |> replace_na(0)
  )
  ) |>
  as.data.frame()

#check that, on a trait-by-trait basis,
# all psnps are more distant than the LD threshold of 270 kb
map(gwas_blink_list_psnp, \(x) { x |>
    arrange(Position_G) |>
    filter(lead(Position_G) - Position_G < 2.7e5)
}
)
#yes, so no need to process further

upset(gwas_blink_psnp_snps,
      intersect = c("BLINK_density", "BLINK_length", "BLINK_width",
                    "BLINK_cvars", "BLINK_field", "BLINK_gh", "BLINK_meta"
      )
)
#essentially consistent with snps

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
    ) |>
    arrange(Chromosome, Position)
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

#add meta-analysis outcomes
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
#generate absolute positions on genome
Chr_max_mlm <- gwas_mlm_list[[1]] |>
  arrange(Chromosome, Position) |>
  group_by(Chromosome) |>
  mutate(Chromosome_max = max(Position)) |>
  ungroup() |>
  distinct(Chromosome, Chromosome_max) |>
  mutate(Chromosome_max_Genome = cumsum(Chromosome_max))

gwas_mlm_list_mta <- gwas_mlm_list |>
  map(\(x) full_join(x, Chr_max_mlm) |>
        arrange(Chromosome, Position) |>
        mutate(Position_G = Chromosome_max_Genome - (Chromosome_max - Position),
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
  xlim(7.67e7, 7.83e7)

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

#check lengths
map_vec(gwas_mlm_list_psnp, nrow)

#check that, on a trait-by-trait basis,
# all psnps are more distant than the LD threshold of 270 kb
map(gwas_mlm_list_psnp, \(x) { x |>
    arrange(Position_G) |>
    filter(lead(Position_G) - Position_G < 2.7e5)
}
)
#no, many of these are in linkage, so I need to filter them out

#function to identify which are true, LD corrected P-SNPs
keep_psnps <- \(x){
  if(nrow(x) == 1){ return(x) } else {
    x <- mutate(x,
                psnp = NA,
                Position_G_plus270 = Position_G + 2.7e5,
                Position_G_minus270 = Position_G - 2.7e5
    )
    print(x)
    #print(any(x$psnp == "psnp", na.rm = TRUE))
    repeat {
      if (all(is.na(x$psnp))) {
        PGnext <- filter(x, neg_log10p == max(neg_log10p)) |>
          select(Position_G) |>
          slice_head() |> #needed to prevent multiples with same neg_log10p
          as.vector()
      } else {
        PGnext <- filter(x, is.na(psnp)) |>
          filter(neg_log10p == max(neg_log10p)) |>
          select(Position_G) |>
          slice_head() |> #needed to prevent multiples with same neg_log10p
          as.vector()
      }
      
      print(PGnext)
      x <- mutate(x,
                  psnp = case_when(
                    psnp == "psnp" | Position_G == PGnext ~ "psnp",
                    .default = NA)
      )
      
      print(x)
      
      #realised this was unnecessary, because it just finds PGnext
      # PG <- filter(x, psnp == "psnp") |>
      #   filter(neg_log10p == min(neg_log10p)) |>
      #   select(Position_G) |>
      #   as.vector()
      
      x <- filter(x,
                  !is.na(psnp) |
                    (Position_G_plus270 < PGnext &
                       Position_G_minus270 < PGnext) |
                    (Position_G_plus270 > PGnext &
                       Position_G_minus270 > PGnext)
      )
      
      if (all(!is.na(x$psnp))) { break }
    }
    
    return(x)
    
  }
  
}

gwas_mlm_list_psnp2 <- gwas_mlm_list_psnp |>
  map_depth(1, \(x) split(x, ~ Chromosome)) |>
  map_depth(2, keep_psnps) |>
  map_depth(1, \(x) list_rbind(x, names_to = "Chromosome"))

map_vec(gwas_mlm_list_psnp, nrow)
map_vec(gwas_mlm_list_psnp2, nrow)

#reformat and use complexupset
gwas_mlm_sig_snps <- gwas_mlm_list_sig |>
  list_rbind(names_to = "method_phenotype") |>
  select(-c(P.value:neg_log10p)) |>
  pivot_wider(names_from = method_phenotype, values_from = mta) |>
  mutate(across(starts_with("MLM"),
                ~ ifelse(.x == "sig", 1, 0) |> replace_na(0)
  )
  )

head(gwas_mlm_sig_snps)
#need to manually reconstruct empty classes
gwas_mlm_sig_snps[ , c("MLM_length", "MLM_width", "MLM_field")] <- 0

upset(gwas_mlm_sig_snps,
      intersect = c("MLM_density", "MLM_length", "MLM_width",
                    "MLM_cvars", "MLM_field", "MLM_gh", "MLM_meta"
      )
)

gwas_mlm_psnp_snps <- gwas_mlm_list_psnp2 |>
  list_rbind(names_to = "method_phenotype") |>
  select(-c(P.value:neg_log10p, mta:max_neg_log10p)) |>
  pivot_wider(names_from = method_phenotype, values_from = psnp) |>
  mutate(across(starts_with("MLM"),
                ~ ifelse(.x == "psnp", 1, 0) |> replace_na(0)
  )
  )

head(gwas_mlm_psnp_snps)
#need to manually reconstruct empty classes
gwas_mlm_psnp_snps[ , c("MLM_length", "MLM_width", "MLM_field")] <- 0

upset(gwas_mlm_psnp_snps,
      intersect = c("MLM_density", "MLM_length", "MLM_width",
                    "MLM_cvars", "MLM_field", "MLM_gh", "MLM_meta"
      )
)

mlm_bnk_sig <- full_join(
  gwas_mlm_sig_snps |>
    pivot_longer(starts_with("MLM"),
                 names_to = "method_phenotype",
                 values_to = "yes_no"
    ),
  gwas_blink_sig_snps |>
    pivot_longer(starts_with("BLINK"),
                 names_to = "method_phenotype",
                 values_to = "yes_no"
    )
)

#write a file with all the mtas listed
write.csv(mlm_bnk_sig, "mlm_blink_mta.csv", row.names = FALSE)

mlm_bnk_sig_wide <- mlm_bnk_sig |>
  pivot_wider(names_from = method_phenotype, values_from = yes_no) |>
  mutate(across(matches("^MLM|^BLINK"), ~ replace_na(.x, 0)))


#To enable pretty plotting, remove underscores in names
names(mlm_bnk_sig_wide) <- str_replace_all(names(mlm_bnk_sig_wide), "_", " ")

a <- upset(mlm_bnk_sig_wide,
           intersect = c("BLINK density", "BLINK length", "BLINK width",
                         "BLINK cvars", "BLINK field", "BLINK gh", "BLINK meta",
                         "MLM density", "MLM length", "MLM width",
                         "MLM cvars", "MLM field", "MLM gh", "MLM meta"
           ),
           sort_sets = FALSE,
           sort_intersections = "descending",
           sort_intersections_by = "cardinality",
           base_annotations = list(
             'Intersection size'= (
               intersection_size()
               + ylim(c(0, 135))
             )
           ),
           set_sizes = upset_set_size() + ylim(170, 0)
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
mlm_bnk_psnp <- full_join(
  gwas_mlm_psnp_snps |>
    select(-contains("270")) |>
    pivot_longer(starts_with("MLM"),
                 names_to = "method_phenotype",
                 values_to = "yes_no"
    ),
  gwas_blink_psnp_snps |>
    pivot_longer(starts_with("BLINK"),
                 names_to = "method_phenotype",
                 values_to = "yes_no"
    )
) |>
  pivot_wider(names_from = method_phenotype, values_from = yes_no) |>
  mutate(across(matches("^MLM|^BLINK"), ~ replace_na(.x, 0)))  

distinct(mlm_bnk_psnp, SNP, Chromosome, Position_G) |>
  arrange(Chromosome, Position_G)

names(mlm_bnk_psnp) <- str_replace_all(names(mlm_bnk_psnp), "_", " ")

b <- upset(mlm_bnk_psnp,
           intersect = c("BLINK density", "BLINK length", "BLINK width",
                         "BLINK cvars", "BLINK field", "BLINK gh", "BLINK meta",
                         "MLM density", "MLM length", "MLM width",
                         "MLM cvars", "MLM field", "MLM gh", "MLM meta"
           ),
           sort_sets = FALSE,
           sort_intersections = "descending",
           sort_intersections_by = "cardinality",
           base_annotations = list(
             'Intersection size'= (
               intersection_size()
               + ylim(c(0, 13.5))
             )
           ),
           set_sizes = upset_set_size() + ylim(170, 0)
)

b

pdf("mlm_blink_upset.pdf", w = 320/25.8, h = 320/25.8)

a / b

dev.off()

mlm_bnk_psnp_table <- list_rbind(
  list(
    list_rbind(gwas_blink_list_psnp, names_to = "method_phenotype"),
    list_rbind(gwas_mlm_list_psnp2, names_to = "method_phenotype") |>
      select(-contains("270"))
  )
)

write.csv(mlm_bnk_psnp_table, "mlm_blink_psnp.csv", row.names = FALSE)

