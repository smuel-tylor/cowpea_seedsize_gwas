#function that filters a set of significant snps to identify those with the
# maximum neg_log10p in any given 270 kb window (psnps).
# Starts by finding and flagging the psnp with the highest neg_log10p
# eliminates all snps within +/- 270 kb of the psnp
# then searches for the next psnp
keep_psnps <- \(x){
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
        as.vector()
    } else {
      PGnext <- filter(x, is.na(psnp)) |>
        filter(neg_log10p == max(neg_log10p)) |>
        select(Position_G) |>
        as.vector()
    }
    
    print(PGnext)
    x <- mutate(x, psnp = case_when(
      psnp == "psnp" | Position_G == PGnext ~ "psnp",
      .default = NA)
    )
    
    print(x)
    
    PG <- filter(x, psnp == "psnp") |>
      filter(neg_log10p == min(neg_log10p)) |>
      select(Position_G) |>
      as.vector()
    
    x <- filter(x,
                !is.na(psnp) |
                  (Position_G_plus270 < PG & Position_G_minus270 < PG) |
                  (Position_G_plus270 > PG & Position_G_minus270 > PG)
    )
    
    if (all(!is.na(x$psnp))) { break }
  }
  
  return(x)
  
}

gwas_blink_list_psnp <- gwas_blink_list_mta |>
  map(\(x) filter(x, mta == "sig")) |>
  map(keep_psnps)

map_vec(gwas_blink_list_mta |>
          map(\(x) filter(x, mta == "sig")), nrow)
map_vec(gwas_blink_list_psnp, nrow)