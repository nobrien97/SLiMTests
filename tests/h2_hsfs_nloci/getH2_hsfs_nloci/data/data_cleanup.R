[nb9894@gadi-cpu-clx-0680 h2_hsfs]$ R    

R version 4.0.0 (2020-04-24) -- "Arbor Day"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> library(tidyverse)
-- Attaching packages ------------------------------------------------------------------------------------------------------ tidyverse 1.3.1 --
v ggplot2 3.3.6     v purrr   0.3.4
v tibble  3.1.7     v dplyr   1.0.9
v tidyr   1.2.0     v stringr 1.4.0
v readr   2.1.2     v forcats 0.5.1
-- Conflicts --------------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
x dplyr::filter() masks stats::filter()
x dplyr::lag()    masks stats::lag()
> d_
> d_qg <- read_csv("slim_qg.csv", col_names = F)
Rows: 201600 Columns: 11                                                                                           
-- Column specification -----------------------------------------------------------------------------------------------------------------------
Delimiter: ","
dbl (11): X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11

i Use `spec()` to retrieve the full column specification for this data.
i Specify the column types or set `show_col_types = FALSE` to quiet this message.
> d_h2 <- read_csv("./getH2_hsfs/out_h2.csv", col_names = F)
Rows: 21959 Columns: 18                                                                                            
-- Column specification -----------------------------------------------------------------------------------------------------------------------
Delimiter: ","
dbl (18): X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, ...

i Use `spec()` to retrieve the full column specification for this data.
i Specify the column types or set `show_col_types = FALSE` to quiet this message.
> names(d_h2) <- c("gen", "seed", "modelindex", "VarA", "VarD", "VarAA", "VarR", 
         "Va+                  "VarA.SE", "VarD.SE", "VarAA.SE", "VarR.SE", "H2.A.Estimate", 
+                  "H2.A.SE", "H2.D.Estimate", "H2.D.SE", "H2.AA.Estimate", 
+                  "H2.AA.SE", "AIC")
> names(d_qg) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", 
+                  "phenovar", "dist", "w", "deltaPheno", "deltaw")
> d_muts <- fre
freeny     freeny.x   freeny.y   frequency  
> d_muts <- fread("d_muts.csv", colClasses = c(seed = 'character'), fill = T)
Error in fread("d_muts.csv", colClasses = c(seed = "character"), fill = T) : 
  could not find function "fread"
> d_muts <- data.table::fread("slim_muts.csv", colClasses = c(seed = 'character'), fill = T)
|--------------------------------------------------|
|==================================================|
Warning message:
In data.table::fread("slim_muts.csv", colClasses = c(seed = "character"),  :
  Column name 'seed' (colClasses[[1]][1]) not found
> d_muts <- read_csv("slim_muts.csv", col_names = F)                                        
                                                                                                                   
Rows: 72645742 Columns: 13
-- Column specification -----------------------------------------------------------------------------------------------------------------------
Delimiter: ","
dbl (12): X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12
lgl  (1): X13

i Use `spec()` to retrieve the full column specification for this data.
i Specify the column types or set `show_col_types = FALSE` to quiet this message.
> 
> head(d_muts)
# A tibble: 6 x 13
     X1         X2    X3    X4    X5    X6    X7    X8     X9   X10    X11   X12
  <dbl>      <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl>
1   500 1002521450     0     3     3    99     0     2 0.849    Inf 0.0014    14
2   500 1002521450     0     3   197   148     0     4 0.221    Inf 0.0573   573
3   500 1002521450     0     3   459   297     0     7 0.269    Inf 0.042    420
4   500 1002521450     0     3   447   386     0     7 0.368    Inf 0.0456   456
5   500 1002521450     0     3   968   233     0    13 0.257    Inf 0.0576   576
6   500 1002521450     0     3   970   505     0    13 0.0299   Inf 0.0366   366
# ... with 1 more variable: X13 <lgl>
> head(d_muts$X13)
[1] NA NA NA NA NA NA
> length(d_muts$X13[!is.na(d_muts$X13)])
[1] 0
> names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "chi", "Freq", "mutCount", "fixGen")
> head(d_qg)
# A tibble: 6 x 11
    gen       seed modelindex   meanH    VA phenomean phenovar  dist     w
  <dbl>      <dbl>      <dbl>   <dbl> <dbl>     <dbl>    <dbl> <dbl> <dbl>
1   500 1002521450          0 0.00811 0.426     0.867    0.426 0.523 0.979
2  1000 1002521450          0 0.0144  0.413     0.922    0.413 0.504 0.980
3  1500 1002521450          0 0.0180  0.344     0.931    0.344 0.444 0.983
4  2000 1002521450          0 0.0214  0.277     0.952    0.277 0.379 0.987
5  2500 1002521450          0 0.0235  0.217     0.955    0.217 0.337 0.989
6  3000 1002521450          0 0.0285  0.216     1.03     0.216 0.325 0.990
# ... with 2 more variables: deltaPheno <dbl>, deltaw <dbl>
> d_muts$seed <- as.factor(d_muts$seed)
index <- as.factor(d_muts$modelindex)> d_muts$modelindex <- as.factor(d_muts$modelindex)
> d_qg$modelindex <- as.factor(d_qg$modelindex)
> d_qg$seed <- as.factor(d_qg$seed)
> d_combined <- inner_join(d_muts, d_qg, by = c("gen", "seed", "modelindex"))
> d_h2$seed <- as.factor(d_h2$seed)
> d_h2$modelindex <- as.factor(d_h2$modelindex)
> d_combined_full <- full_join(d_combined, d_h2, by = c("gen", "seed", "modelindex"))
> head(d_combined_full) 
# A tibble: 6 x 36
    gen seed       modelindex mutType mutID position constraint originGen  value
  <dbl> <fct>      <fct>        <dbl> <dbl>    <dbl>      <dbl>     <dbl>  <dbl>
1   500 1002521450 0                3     3       99          0         2 0.849 
2   500 1002521450 0                3   197      148          0         4 0.221 
3   500 1002521450 0                3   459      297          0         7 0.269 
4   500 1002521450 0                3   447      386          0         7 0.368 
5   500 1002521450 0                3   968      233          0        13 0.257 
6   500 1002521450 0                3   970      505          0        13 0.0299
# ... with 27 more variables: chi <dbl>, Freq <dbl>, mutCount <dbl>,
#   fixGen <lgl>, meanH <dbl>, VA <dbl>, phenomean <dbl>, phenovar <dbl>,
#   dist <dbl>, w <dbl>, deltaPheno <dbl>, deltaw <dbl>, VarA <dbl>,
#   VarD <dbl>, VarAA <dbl>, VarR <dbl>, VarA.SE <dbl>, VarD.SE <dbl>,
#   VarAA.SE <dbl>, VarR.SE <dbl>, H2.A.Estimate <dbl>, H2.A.SE <dbl>,
#   H2.D.Estimate <dbl>, H2.D.SE <dbl>, H2.AA.Estimate <dbl>, H2.AA.SE <dbl>,
#   AIC <dbl>
> d_combined_full[1,]  
# A tibble: 1 x 36
    gen seed  modelindex mutType mutID position constraint originGen value   chi
  <dbl> <fct> <fct>        <dbl> <dbl>    <dbl>      <dbl>     <dbl> <dbl> <dbl>
1   500 1002~ 0                3     3       99          0         2 0.849   Inf
# ... with 26 more variables: Freq <dbl>, mutCount <dbl>, fixGen <lgl>,
#   meanH <dbl>, VA <dbl>, phenomean <dbl>, phenovar <dbl>, dist <dbl>,
#   w <dbl>, deltaPheno <dbl>, deltaw <dbl>, VarA <dbl>, VarD <dbl>,
#   VarAA <dbl>, VarR <dbl>, VarA.SE <dbl>, VarD.SE <dbl>, VarAA.SE <dbl>,
#   VarR.SE <dbl>, H2.A.Estimate <dbl>, H2.A.SE <dbl>, H2.D.Estimate <dbl>,
#   H2.D.SE <dbl>, H2.AA.Estimate <dbl>, H2.AA.SE <dbl>, AIC <dbl>
> d_combined_full[1,] %>% print(n=36)
# A tibble: 1 x 36
    gen seed  modelindex mutType mutID position constraint originGen value   chi
  <dbl> <fct> <fct>        <dbl> <dbl>    <dbl>      <dbl>     <dbl> <dbl> <dbl>
1   500 1002~ 0                3     3       99          0         2 0.849   Inf
# ... with 26 more variables: Freq <dbl>, mutCount <dbl>, fixGen <lgl>,
#   meanH <dbl>, VA <dbl>, phenomean <dbl>, phenovar <dbl>, dist <dbl>,
#   w <dbl>, deltaPheno <dbl>, deltaw <dbl>, VarA <dbl>, VarD <dbl>,
#   VarAA <dbl>, VarR <dbl>, VarA.SE <dbl>, VarD.SE <dbl>, VarAA.SE <dbl>,
#   VarR.SE <dbl>, H2.A.Estimate <dbl>, H2.A.SE <dbl>, H2.D.Estimate <dbl>,
#   H2.D.SE <dbl>, H2.AA.Estimate <dbl>, H2.AA.SE <dbl>, AIC <dbl>
> d_combined_full[1,] %>% as_tibble() %>% print(n=36)
# A tibble: 1 x 36
    gen seed  modelindex mutType mutID position constraint originGen value   chi
  <dbl> <fct> <fct>        <dbl> <dbl>    <dbl>      <dbl>     <dbl> <dbl> <dbl>
1   500 1002~ 0                3     3       99          0         2 0.849   Inf
# ... with 26 more variables: Freq <dbl>, mutCount <dbl>, fixGen <lgl>,
#   meanH <dbl>, VA <dbl>, phenomean <dbl>, phenovar <dbl>, dist <dbl>,
#   w <dbl>, deltaPheno <dbl>, deltaw <dbl>, VarA <dbl>, VarD <dbl>,
#   VarAA <dbl>, VarR <dbl>, VarA.SE <dbl>, VarD.SE <dbl>, VarAA.SE <dbl>,
#   VarR.SE <dbl>, H2.A.Estimate <dbl>, H2.A.SE <dbl>, H2.D.Estimate <dbl>,
#   H2.D.SE <dbl>, H2.AA.Estimate <dbl>, H2.AA.SE <dbl>, AIC <dbl>
> d_combined_full[1,] %>% as_tibble() %>% print(width=36)
# A tibble: 1 x 36
    gen seed      modelindex mutType
  <dbl> <fct>     <fct>        <dbl>
1   500 10025214~ 0                3
# ... with 32 more variables:
#   mutID <dbl>, position <dbl>,
#   constraint <dbl>,
#   originGen <dbl>, value <dbl>,
#   chi <dbl>, Freq <dbl>,
#   mutCount <dbl>, fixGen <lgl>,
#   meanH <dbl>, VA <dbl>, ...
> d_combined_full[1,] %>% as_tibble() %>% print(width=Inf)
# A tibble: 1 x 36
    gen seed       modelindex mutType mutID position constraint originGen value
  <dbl> <fct>      <fct>        <dbl> <dbl>    <dbl>      <dbl>     <dbl> <dbl>
1   500 1002521450 0                3     3       99          0         2 0.849
    chi   Freq mutCount fixGen   meanH    VA phenomean phenovar  dist     w
  <dbl>  <dbl>    <dbl> <lgl>    <dbl> <dbl>     <dbl>    <dbl> <dbl> <dbl>
1   Inf 0.0014       14 NA     0.00811 0.426     0.867    0.426 0.523 0.979
  deltaPheno deltaw  VarA  VarD VarAA  VarR VarA.SE VarD.SE VarAA.SE VarR.SE
       <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
1      0.523  0.979    NA    NA    NA    NA      NA      NA       NA      NA
  H2.A.Estimate H2.A.SE H2.D.Estimate H2.D.SE H2.AA.Estimate H2.AA.SE   AIC
          <dbl>   <dbl>         <dbl>   <dbl>          <dbl>    <dbl> <dbl>
1            NA      NA            NA      NA             NA       NA    NA
> d_combined_full[800000,] %>% as_tibble() %>% print(width=Inf)
# A tibble: 1 x 36
    gen seed       modelindex mutType mutID position constraint originGen value
  <dbl> <fct>      <fct>        <dbl> <dbl>    <dbl>      <dbl>     <dbl> <dbl>
1 51350 1070070213 1                3  2002      901          0        24 0.549
            chi  Freq mutCount fixGen meanH    VA phenomean phenovar    dist
          <dbl> <dbl>    <dbl> <lgl>  <dbl> <dbl>     <dbl>    <dbl>   <dbl>
1 0.00000000862     1    10000 NA     0.115  46.5   2.03e15  2.05e34 2.03e15
      w deltaPheno    deltaw     VarA  VarD VarAA    VarR VarA.SE VarD.SE
  <dbl>      <dbl>     <dbl>    <dbl> <dbl> <dbl>   <dbl>   <dbl>   <dbl>
1 0.816    2.03e15 -0.000436 -4.05e14     0     0 6.98e16 7.54e14 1.65e14
  VarAA.SE VarR.SE H2.A.Estimate H2.A.SE H2.D.Estimate H2.D.SE H2.AA.Estimate
     <dbl>   <dbl>         <dbl>   <dbl>         <dbl>   <dbl>          <dbl>
1  2.80e15 4.30e15      -0.00584  0.0109             0 0.00238              0
  H2.AA.SE   AIC
     <dbl> <dbl>
1   0.0403 1000.
> data.table::fwrite(d_combined_full, "d_combined_full.csv")
> d_combined_after <- d_combined_full %>% filter(gen >= 50000)                                                                       
> nrow(d_combined_after)
[1] 26242585
> nrow(d_combined_full)
[1] 73197498
> data.table::fwrite(d_combined_after, "d_combined_after.csv")