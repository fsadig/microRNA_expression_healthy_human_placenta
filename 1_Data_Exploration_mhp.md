---
title: "SSRP: Data Cleanup + Exploration"
author: Nikita Telkar
output: 
  html_document:
    keep_md: yes
  # rmdformats::downcute:
  #   self_contained: true
  #   default_style: "light"
  #   downcute_theme: "default"
    highlight: haddock #tango, pygments, kate, monochrome, espresso, zenburn, haddock, textmate
    toc: true 
    toc_depth: 4
    toc_float: 
      collapsed: true 
      smooth_scroll: false  
    theme: sandstone  #cosmo, paper, lumen, sandstone, simplex, yeti, cerulean, journal, flatly, darkly, readable, spacelab, united

---


***  

### 0.0 Introduction  

Here, we're going to go through what our raw data and clinical data looks like, and then clean or 'wrangle' our data to a useable format. We're also going to make some plots depicting the spread of our data.   

***  

### 1.0 Loading Packages and Data


```r
# rmarkdown packages
library(knitr) 
library(rmarkdown)
library(devtools)
```

```
## Loading required package: usethis
```

```r
# data wrangling packages
library(tidyverse)
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
```

```
## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
## ✔ tibble  3.1.6     ✔ dplyr   1.0.8
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(reshape2)
```

```
## 
## Attaching package: 'reshape2'
```

```
## The following object is masked from 'package:tidyr':
## 
##     smiths
```

```r
library(magrittr)
```

```
## 
## Attaching package: 'magrittr'
```

```
## The following object is masked from 'package:purrr':
## 
##     set_names
```

```
## The following object is masked from 'package:tidyr':
## 
##     extract
```

```r
# data loading packages
library(readxl) 
library(openxlsx)

# plotting/table packages
library(arsenal)
```

```
## 
## Attaching package: 'arsenal'
```

```
## The following object is masked from 'package:magrittr':
## 
##     set_attr
```

```r
library(kableExtra) 
```

```
## 
## Attaching package: 'kableExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     group_rows
```

```r
library(egg)
```

```
## Loading required package: gridExtra
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```r
library(RColorBrewer)
library(ggpubr)
```

```
## 
## Attaching package: 'ggpubr'
```

```
## The following object is masked from 'package:egg':
## 
##     ggarrange
```

```r
library(emo)

#BiocManager::install(c("limma", "edgeR", "DESeq2", "biomaRt", "miRBaseConverter", "multiMiR", "preprocessCore"))
#remotes::install_github("wvictor14/plomics")
```

We're going to use a very handy package called the `here` package to navigate between our main folder, and to read in and save files  


```r
library(here)
```

```
## here() starts at /Users/fidan/BCCHR/SSRP
```

The output should show you which folder the `here` package sets as its root/main directory.  

Now, read in the 2 excel files and the 1 miRNA data file located within your `data` folder


```r
lam_ids <- read_excel(here::here("data", "LamLab_IDs.xlsx"))
rob_ids <- read_excel (here:: here ("data", "RobLab_IDs.xlsx"))
```

```
## New names:
## * `` -> ...1
```

```r
mirna <- read_excel(here:: here("data", "miRNA_quantification_raw.xlsx"))
```

Can you find the function we can use to view the dimensions, meaning the number of rows and columns, of the dataframe? How about to view the structure of the dataframe, meaning the types of variables (e.g., words/characters, numeric, etc)?  

*Hint: there's two functions to view the structure, one is a baseR function, and one is from the tidyverse package.*    


```r
dim(lam_ids)
```

```
## [1] 30  4
```

```r
dim(rob_ids)
```

```
## [1] 30 20
```

```r
dim(mirna)
```

```
## [1] 2887   32
```

```r
str(lam_ids)
```

```
## tibble [30 × 4] (S3: tbl_df/tbl/data.frame)
##  $ GSC_ID                 : chr [1:30] "MX1303-C5JC4ACXX-3-TAGGAT" "MX1304-C5JC4ACXX-4-CCGGTG" "MX1305-C5JC4ACXX-5-TGTTGG" "MX1306-C5JC4ACXX-6-GTATAG" ...
##  $ ID                     : chr [1:30] "FT42_cv" "FT52_cv" "FT54_cv" "FT60_cv" ...
##  $ lamID                  : chr [1:30] "Rob42_cv" "Rob52_cv" "Rob54_cv" "Rob60_cv" ...
##  $ totalReadsAligned_miRNA: chr [1:30] "7136633" "10118247" "19024542" "13461776" ...
```

```r
str(rob_ids) 
```

```
## tibble [30 × 20] (S3: tbl_df/tbl/data.frame)
##  $ ...1           : chr [1:30] "1" "2" "3" "4" ...
##  $ ID             : chr [1:30] "FT42_cv" "FT52_cv" "FT54_cv" "FT60_cv" ...
##  $ Case ID        : chr [1:30] "FT42" "FT52" "FT54" "FT60" ...
##  $ Tissue         : chr [1:30] "chorionic villi" "chorionic villi" "chorionic villi" "chorionic villi" ...
##  $ Condition      : chr [1:30] "con" "NTD" "NTD" "NTD" ...
##  $ Extra Condition: chr [1:30] "pPROM" "anencephaly" "spina bifida" "spina bifida" ...
##  $ Gestational Age: chr [1:30] "17" "21.8" "20.399999999999999" "23.4" ...
##  $ Trimester      : chr [1:30] "2" "2" "2" "2" ...
##  $ Sex            : chr [1:30] "MALE" "FEMALE" "MALE" "FEMALE" ...
##  $ Processing_time: chr [1:30] "72" "24" "48" "48" ...
##  $ libraryID      : chr [1:30] "MX1303" "MX1304" "MX1305" "MX1306" ...
##  $ flowCell       : chr [1:30] "C5JC4ACXX" "C5JC4ACXX" "C5JC4ACXX" "C5JC4ACXX" ...
##  $ Lane           : chr [1:30] "3" "4" "5" "6" ...
##  $ indexAdapter   : chr [1:30] "TAGGAT" "CCGGTG" "TGTTGG" "GTATAG" ...
##  $ MLPA           : chr [1:30] "46,XY" "46,XX" "46,XY" "46,XX" ...
##  $ misoprostol    : chr [1:30] "UNK" "UNK" "UNK" "Y" ...
##  $ modeLabour     : chr [1:30] "D+C" "D+E" "IOL" "IOL" ...
##  $ RQS            : chr [1:30] "3.3" "4" "3.6" "3.4" ...
##  $ Degraded       : chr [1:30] "NA" "NA" "NA" "NA" ...
##  $ 450k.array     : chr [1:30] "Match" "Match" "Match" "Match" ...
```

```r
str(mirna)
```

```
## tibble [2,887 × 32] (S3: tbl_df/tbl/data.frame)
##  $ precursor                : chr [1:2887] "hsa-let-7a-1" "hsa-let-7a-1" "hsa-let-7a-2" "hsa-let-7a-2" ...
##  $ miRNA                    : chr [1:2887] "hsa-let-7a-3p" "hsa-let-7a-5p" "hsa-let-7a-2-3p" "hsa-let-7a-5p" ...
##  $ MX1303-C5JC4ACXX-3-TAGGAT: num [1:2887] 444 96079 3 96112 450 ...
##  $ MX1304-C5JC4ACXX-4-CCGGTG: num [1:2887] 647 164316 2 164372 649 ...
##  $ MX1305-C5JC4ACXX-5-TGTTGG: num [1:2887] 709 293191 22 293225 712 ...
##  $ MX1306-C5JC4ACXX-6-GTATAG: num [1:2887] 813 226991 12 227042 823 ...
##  $ MX1307-C5JC4ACXX-7-AGCATC: num [1:2887] 495 186434 8 186507 501 ...
##  $ MX1307-C5JC4ACXX-7-CAGGCC: num [1:2887] 314 55365 3 55392 314 ...
##  $ MX1310-C5JC1ACXX-4-CTCTAC: num [1:2887] 276 63794 0 63813 278 ...
##  $ MX1310-C5JC1ACXX-4-GGAACT: num [1:2887] 1493 522367 16 522475 1506 ...
##  $ MX1310-C5JC1ACXX-4-GGACGG: num [1:2887] 1183 396173 3 396201 1196 ...
##  $ MX1310-C5JC1ACXX-4-TGACAT: num [1:2887] 1483 306320 7 306488 1498 ...
##  $ MX1355-C6RGTANXX-2-AAGCTA: num [1:2887] 37 4892 2 4885 37 ...
##  $ MX1355-C6RGTANXX-2-ACATCG: num [1:2887] 377 34289 2 34249 380 ...
##  $ MX1355-C6RGTANXX-2-CAAGTT: num [1:2887] 26 4281 0 4278 26 ...
##  $ MX1355-C6RGTANXX-2-CATTCA: num [1:2887] 82 18877 1 18851 82 ...
##  $ MX1355-C6RGTANXX-2-GGAACT: num [1:2887] 161 54184 1 54159 163 ...
##  $ MX1356-C6RGTANXX-3-ATGGCA: num [1:2887] 205 53078 5 53113 205 ...
##  $ MX1356-C6RGTANXX-3-CCTTGC: num [1:2887] 373 97396 2 97363 374 ...
##  $ MX1356-C6RGTANXX-3-CGGCCT: num [1:2887] 63 8462 1 8458 63 ...
##  $ MX1356-C6RGTANXX-3-GCGTGG: num [1:2887] 284 36696 0 36658 287 ...
##  $ MX1356-C6RGTANXX-3-GCTGTA: num [1:2887] 801 171395 3 171415 803 ...
##  $ MX1356-C6RGTANXX-3-GTATAG: num [1:2887] 717 154473 9 154346 723 ...
##  $ MX1356-C6RGTANXX-3-TAGTTG: num [1:2887] 43 3864 3 3867 44 ...
##  $ MX1356-C6RGTANXX-3-TGACAT: num [1:2887] 1761 382620 5 382570 1765 ...
##  $ MX1357-C6RGTANXX-4-AATTAT: num [1:2887] 146 31237 0 31180 146 ...
##  $ MX1357-C6RGTANXX-4-AGTCTT: num [1:2887] 737 186118 8 185711 745 ...
##  $ MX1357-C6RGTANXX-4-CATGGG: num [1:2887] 120 27350 0 27281 122 ...
##  $ MX1357-C6RGTANXX-4-GCCTAA: num [1:2887] 845 161419 5 161003 851 ...
##  $ MX1357-C6RGTANXX-4-GTAGCC: num [1:2887] 613 109686 13 109434 614 ...
##  $ MX1357-C6RGTANXX-4-TATCGT: num [1:2887] 1103 214843 4 214606 1104 ...
##  $ MX1357-C6RGTANXX-4-TCTGAG: num [1:2887] 944 186632 4 186168 959 ...
```

```r
## not sure about the tidyverse one
```


Let's now look at the first few rows of each of the files, or dataframes, that we've loaded and are going to use. I've written down the function using the `lam_ids` datafram, or df, as an example.  

**Do the same for the other dfs**


```r
head(lam_ids)
```

```
## # A tibble: 6 × 4
##   GSC_ID                    ID      lamID    totalReadsAligned_miRNA
##   <chr>                     <chr>   <chr>    <chr>                  
## 1 MX1303-C5JC4ACXX-3-TAGGAT FT42_cv Rob42_cv 7136633                
## 2 MX1304-C5JC4ACXX-4-CCGGTG FT52_cv Rob52_cv 10118247               
## 3 MX1305-C5JC4ACXX-5-TGTTGG FT54_cv Rob54_cv 19024542               
## 4 MX1306-C5JC4ACXX-6-GTATAG FT60_cv Rob60_cv 13461776               
## 5 MX1307-C5JC4ACXX-7-AGCATC FT67_cv Rob67_cv 15587802               
## 6 MX1307-C5JC4ACXX-7-CAGGCC FT73_cv Rob73_cv 9102519
```

```r
head(rob_ids)
```

```
## # A tibble: 6 × 20
##   ...1  ID      `Case ID` Tissue     Condition `Extra Conditi…` `Gestational A…`
##   <chr> <chr>   <chr>     <chr>      <chr>     <chr>            <chr>           
## 1 1     FT42_cv FT42      chorionic… con       pPROM            17              
## 2 2     FT52_cv FT52      chorionic… NTD       anencephaly      21.8            
## 3 3     FT54_cv FT54      chorionic… NTD       spina bifida     20.399999999999…
## 4 4     FT60_cv FT60      chorionic… NTD       spina bifida     23.4            
## 5 5     FT67_cv FT67      chorionic… NTD       spina bifida     22.4            
## 6 6     FT73_cv FT73      chorionic… con       pPROM            18              
## # … with 13 more variables: Trimester <chr>, Sex <chr>, Processing_time <chr>,
## #   libraryID <chr>, flowCell <chr>, Lane <chr>, indexAdapter <chr>,
## #   MLPA <chr>, misoprostol <chr>, modeLabour <chr>, RQS <chr>, Degraded <chr>,
## #   `450k.array` <chr>
```

```r
head(mirna)
```

```
## # A tibble: 6 × 32
##   precursor    miRNA          `MX1303-C5JC4A…` `MX1304-C5JC4A…` `MX1305-C5JC4A…`
##   <chr>        <chr>                     <dbl>            <dbl>            <dbl>
## 1 hsa-let-7a-1 hsa-let-7a-3p               444              647              709
## 2 hsa-let-7a-1 hsa-let-7a-5p             96079           164316           293191
## 3 hsa-let-7a-2 hsa-let-7a-2-…                3                2               22
## 4 hsa-let-7a-2 hsa-let-7a-5p             96112           164372           293225
## 5 hsa-let-7a-3 hsa-let-7a-3p               450              649              712
## 6 hsa-let-7a-3 hsa-let-7a-5p             96864           165176           294468
## # … with 27 more variables: `MX1306-C5JC4ACXX-6-GTATAG` <dbl>,
## #   `MX1307-C5JC4ACXX-7-AGCATC` <dbl>, `MX1307-C5JC4ACXX-7-CAGGCC` <dbl>,
## #   `MX1310-C5JC1ACXX-4-CTCTAC` <dbl>, `MX1310-C5JC1ACXX-4-GGAACT` <dbl>,
## #   `MX1310-C5JC1ACXX-4-GGACGG` <dbl>, `MX1310-C5JC1ACXX-4-TGACAT` <dbl>,
## #   `MX1355-C6RGTANXX-2-AAGCTA` <dbl>, `MX1355-C6RGTANXX-2-ACATCG` <dbl>,
## #   `MX1355-C6RGTANXX-2-CAAGTT` <dbl>, `MX1355-C6RGTANXX-2-CATTCA` <dbl>,
## #   `MX1355-C6RGTANXX-2-GGAACT` <dbl>, `MX1356-C6RGTANXX-3-ATGGCA` <dbl>, …
```

The `mirna` df shows how expression data usually looks like - your genes are your rows, and your samples are your columns.  

**But here, there's 2 columns for the miRNAs - can you explain why?** 
## The first column is the precursor miRNA which the microRNA before leaving the nucleus and the second column is the mature miRNA which resulted by further manipulations by proteins like Dicer in the cytoplasm.  

***  

### 2.0 Making eDat and pDat  



#### 2.1 pDat

Now, there's two dfs which have some kinds of IDs in them, and it's one ID type that is present as the sample name in the expression df.  

> Having several files containing different aspects of our required data is very common. In most cases, you are the one who will have to go in and clean up the data to fit the format that you prefer

To make one df which contains all of our phenotype information in it, called pDat, we will have to join the two id dfs  

  
But before we join the dfs, R doesn't do well with column/row names which have spaces within the words.  

**So, let's rename the column names which have spaces between the words**    
*Hint: use the `tidyverse::rename` function*. It's best practice to have an underscore `_` as the separator/delimiter between words  


```r
head(rob_ids)
```

```
## # A tibble: 6 × 20
##   ...1  ID      `Case ID` Tissue     Condition `Extra Conditi…` `Gestational A…`
##   <chr> <chr>   <chr>     <chr>      <chr>     <chr>            <chr>           
## 1 1     FT42_cv FT42      chorionic… con       pPROM            17              
## 2 2     FT52_cv FT52      chorionic… NTD       anencephaly      21.8            
## 3 3     FT54_cv FT54      chorionic… NTD       spina bifida     20.399999999999…
## 4 4     FT60_cv FT60      chorionic… NTD       spina bifida     23.4            
## 5 5     FT67_cv FT67      chorionic… NTD       spina bifida     22.4            
## 6 6     FT73_cv FT73      chorionic… con       pPROM            18              
## # … with 13 more variables: Trimester <chr>, Sex <chr>, Processing_time <chr>,
## #   libraryID <chr>, flowCell <chr>, Lane <chr>, indexAdapter <chr>,
## #   MLPA <chr>, misoprostol <chr>, modeLabour <chr>, RQS <chr>, Degraded <chr>,
## #   `450k.array` <chr>
```

```r
rob_ids <- rob_ids %>% 
   rename_with(~ gsub(" ","_", .x), contains(" "))   

#I found this function when I was trying to find a way to make this work rename (Case_ID = Case ID); this wouldn't run as there is space between the Case and ID
```


Using a function from the tidyverse package:  

**join `lamlab_ids` and `roblab_ids` to create a new dataframe called `pDat`**  

You'll also observe that there's an extra serial number column in the roblab_ids df:  

**remove it too.**  

Further:

**rearrange the order of the columns** to make Case_ID as the first column, followed by ID, lamID, GSC_ID, and then the rest of the columns  
*Hint: specify in your code that you want to use the `select` function from the `dplyr` package*  


```r
pDat <- cbind(lam_ids, rob_ids) %>% 
  subset(select = - ...1) %>% 
  dplyr:: select(Case_ID, ID, lamID, GSC_ID, everything()) 
```

Now, check the variable types of each of our variables in pDat  


```r
str(pDat)
```

```
## 'data.frame':	30 obs. of  23 variables:
##  $ Case_ID                : chr  "FT42" "FT52" "FT54" "FT60" ...
##  $ ID                     : chr  "FT42_cv" "FT52_cv" "FT54_cv" "FT60_cv" ...
##  $ lamID                  : chr  "Rob42_cv" "Rob52_cv" "Rob54_cv" "Rob60_cv" ...
##  $ GSC_ID                 : chr  "MX1303-C5JC4ACXX-3-TAGGAT" "MX1304-C5JC4ACXX-4-CCGGTG" "MX1305-C5JC4ACXX-5-TGTTGG" "MX1306-C5JC4ACXX-6-GTATAG" ...
##  $ totalReadsAligned_miRNA: chr  "7136633" "10118247" "19024542" "13461776" ...
##  $ ID.1                   : chr  "FT42_cv" "FT52_cv" "FT54_cv" "FT60_cv" ...
##  $ Tissue                 : chr  "chorionic villi" "chorionic villi" "chorionic villi" "chorionic villi" ...
##  $ Condition              : chr  "con" "NTD" "NTD" "NTD" ...
##  $ Extra_Condition        : chr  "pPROM" "anencephaly" "spina bifida" "spina bifida" ...
##  $ Gestational_Age        : chr  "17" "21.8" "20.399999999999999" "23.4" ...
##  $ Trimester              : chr  "2" "2" "2" "2" ...
##  $ Sex                    : chr  "MALE" "FEMALE" "MALE" "FEMALE" ...
##  $ Processing_time        : chr  "72" "24" "48" "48" ...
##  $ libraryID              : chr  "MX1303" "MX1304" "MX1305" "MX1306" ...
##  $ flowCell               : chr  "C5JC4ACXX" "C5JC4ACXX" "C5JC4ACXX" "C5JC4ACXX" ...
##  $ Lane                   : chr  "3" "4" "5" "6" ...
##  $ indexAdapter           : chr  "TAGGAT" "CCGGTG" "TGTTGG" "GTATAG" ...
##  $ MLPA                   : chr  "46,XY" "46,XX" "46,XY" "46,XX" ...
##  $ misoprostol            : chr  "UNK" "UNK" "UNK" "Y" ...
##  $ modeLabour             : chr  "D+C" "D+E" "IOL" "IOL" ...
##  $ RQS                    : chr  "3.3" "4" "3.6" "3.4" ...
##  $ Degraded               : chr  "NA" "NA" "NA" "NA" ...
##  $ 450k.array             : chr  "Match" "Match" "Match" "Match" ...
```

We see that all of the variables that should have be numeric are actually encoded as characters - Gestational_Age, Processing_time, Lane, RQS, and totalReadsAligned_miRNA. This will give us problems later on.  

**Convert the above variables from character to numeric**  
*Hint: Use the `mutate_at` function to chnage specific variables*  



```r
pDat <- pDat %>% 
  mutate_at (c("Lane", "Gestational_Age","Processing_time","RQS","totalReadsAligned_miRNA"), as.numeric)
```

```
## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion
```

For the sample which got an `NA` value introduced for its RQS,  

**Replace the NA with 0**  


```r
pDat <- pDat %>% 
  mutate_at(21, ~replace_na(.,0))
```




***  

#### 2.2 eDat

Okay, now if you check our `mirna` df, we can see that the GSC_ID names (the names given to the samples by the GSC (the Genome Sciences Centre) which is the sequencing facility where we've sequenced these samples) are actually not intuitive at all.  

But, now that we've matched the id dataframes, we know exactly which of the GSC_IDs correspond to the Case ID variable.  

**Change the GSC_ID column names to Case_ID names**  

We first have to check whether the GSC_ID order of the column names matches the order of the GSC_ID row names in pDat   


```r
colnames(mirna)[3:32] == pDat$GSC_ID
```

```
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```

```r
# if true for all, proceed. If false, do:

# match(lam_ids$GSC_ID, colnames(mirna)[3:32])
# reorder_idx <- match(pDat$GSC_ID, colnames(mirna)[3:32])
# mirna[3:32] <- mirna[reorder_idx]
# colnames(mirna)[3:32] == pDat$GSC_ID
```

Now, renames the column names in our `mirna` df  
*Hint: use the `colnames` function* to by first selecting the columns to renames, and then the variable you want to rename them by  


```r
colnames(mirna)[3:32] <- pDat$Case_ID
```

Now, you know that there's two columns listing miRNA names in the `mirna` df - a lot of our downstream functions require us to only have numeric values in our df, with the gene names as rownames. But, the `column_to_rownames` function only takes the input of one column, but here we have two  

**Combine the precursor and miRNA columns in `mirna` into one column** using `:` as the separator.  


```r
mirna <- mirna %>% 
  unite("Precursor:miRNA", precursor:miRNA, sep = ":", remove= TRUE)
```

Once combined, 

**Set the new combined column as the rownames for the `mirna` df**  


```r
mirna <- mirna %>%
  column_to_rownames(var="Precursor:miRNA")
```

> This is the format an expression data dataframe should look like 

What is if want to get the average expression for each miRNA across all of the samples?  

**Make a new column showing the average expression across each column per row**  
*Hint: use the `rowSums` function*  


```r
mirna <- mirna %>%
  mutate(Sums = rowSums(.))
```

Let's now make our final expression dataframe named `eDat`  

**Make eDat omitting the last expression sum column from the `mirna` df**  


```r
eDat <- select(mirna, -Sums)
```

***  

### 3.0 Data Exploration  

Okay, now that we have both our `eDat` and `pDat` dfs properly formatted, let's make some plots to display the structure of our dataset  

`ggplot2` (the Grammar of Graphics) package is used to make graphs in R, where it has the following necessary components:  

- `aes` = aesthetics = your x and y axes, and any additional variables you want to deliniate your data by  
- `geoms` = the type to graph you want to plot (e.g., bar plot = geom_bar, scatter plot = geom_point)  

If you look at our eDat df, you'll see that we have each sample as a different column - but if we want to plot our x axis as our samples, ggplot only takes one variable name (column name) for the x axis

And so, we have to first reformat eDat to a df that is plot-able by ggplot. We'll use the slightly advanced function `pivot_longer` to convert our wide eDat df, into a 'longer' df. When this function is used, R forgets to incorporate the rownames into the new df that we make, and hence, we have to change our rownames back into separate columns  

**Convert the rownames of `eDat` into the original two separate columns in a new df**  


```r
eDat_nonames <- eDat %>% 
  rownames_to_column(var="Precursor:miRNA") %>% 
  separate("Precursor:miRNA", c("precursor", "miRNA"), ":")
```

Once we have our original format of the expression dataframe (how it looked when we first loaded it in), we'll now convert it to a longer format  


```r
eDat_plot <- eDat_nonames %>% 
  pivot_longer(cols = -c(precursor, miRNA), names_to = "Case_ID", values_to = "Reads")
```

You'll see that all of the column names (`names_to`) have been gathered up into one column named Case_ID (to match the name of the column with all the samples in pDat), and all the expression counts (`values_to`) have been formatted into one column I've names as Reads  

To add all the phenotype and expression data together in one df to be able to plot by different variables:  

**Join `eDat_plot` with `pDat`**  


```r
eDat_plot <- eDat_plot %>% 
  inner_join( pDat, c ("Case_ID"))
```


Okay, now using pDat, let's make a very simple bar counting the number of samples we have for each 


```r
pDat %>% 
  ggplot(aes(x = Trimester)) +
    #the plus sign is how you add additional arguments/components to your ggplot
  geom_bar(stat = "count") 
```

![](1_Data_Exploration_files/figure-html/tri-plot-1-1.png)<!-- -->

```r
    #here, we're not adding a y axis, because we simply want to count the number of samples, and hence specify "count"
```

Let's add some colours to distinguish our samples, and make the width of the bars smaller  


```r
pDat %>% 
  ggplot(aes(x = Trimester, fill = Trimester)) +
  geom_bar(stat = "count", width = 0.6) 
```

![](1_Data_Exploration_files/figure-html/tri-plot-2-1.png)<!-- -->

Adding labels to the plot, and changing the background:  


```r
pDat %>% 
  ggplot(aes(x = Trimester, fill = Trimester)) +
  geom_bar(stat = "count", width = 0.6) +
  theme_minimal() +
  labs(title = "Distribution of Samples by Trimester", subtitle = "n = 30", x = "Trimester", y = "Number of Samples")
```

![](1_Data_Exploration_files/figure-html/tri-plot-3-1.png)<!-- -->

**Using pDat, I'd like you to plot:**  

- A bar plot showing distribution of samples by sex  


```r
pDat %>% 
  ggplot(aes(x = Sex, fill = Sex)) +
  geom_bar(stat = "count", width = 0.6) +
  theme_minimal() +
  labs(title = "Distribution of Samples by Sex", subtitle = "n = 30", x = "Sex", y = "Number of Samples")
```

![](1_Data_Exploration_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

- A bar plot showing the total reads of aligned miRNAs per sample  


```r
pDat %>% 
  ggplot(aes(x = Case_ID, y = totalReadsAligned_miRNA, fill = Case_ID)) +
  geom_bar (stat = "identity", width = 0.6) +
  theme_minimal(base_size = 4) +
  labs(title = "Distribution of Aligned miRNAs per Sample", x = "Case ID", y = "Number of aligned miRNA")
```

![](1_Data_Exploration_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

- A scatter plot with processing time on the x-axis, and RQS on the y-axis. Is there a correlation between these two variables?  



                    
                    

```r
result <-cor.test(x= pDat$RQS,y= pDat$Processing_time, 
                    method = "pearson")

#There is moderate negative correlation of -0.6027559 with 0.0004233 p value. 

pDat %>% 
  ggplot( aes(x=Processing_time, y=RQS)) +
  geom_point() +
  labs(title = "Processing Time vs RQS", subtitle = result$estimate, x = "Processing Time ", y = "RQS")
```

![](1_Data_Exploration_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
### Couldn't put this line: estimate: result$estimate
```

**Using eDat_plot, I'd like you to plot:**  

by converting Reads to log2(Reads)  

- A scatter plot showing the expression of hsa-miR-100-3p for each sample, colored by trimester *(Hint: use the `filter` function)*  


```r
eDat_plot %>% 
  mutate(Reads=log2(Reads)) %>% 
  filter(miRNA == "hsa-miR-100-3p") %>% 
  ggplot( aes(x=Case_ID, y= Reads, color = Trimester )) +
  geom_point() +
  theme_minimal(base_size = 5) +
  labs(title = "Change of expression of hsa-miR-100-3p in samples over trimesters", x = "Case_ID ", y = "log 2 fold of Reads")
```

![](1_Data_Exploration_files/figure-html/unnamed-chunk-4-1.png)<!-- -->




- A box plot showing the expressions of hsa-miR-21-5p, hsa-miR-200a-5p, and hsa-miR-373-5p, fill by trimester 

```r
eDat_plot %>% 
  mutate(Reads=log2(Reads)) %>% 
  filter (miRNA == c("hsa-miR-21-5p","hsa-miR-200a-5p","hsa-miR-373-5p")) %>% 
  ggplot(aes(x = miRNA, y = Reads, fill = Trimester)) +
  geom_bar (stat = "identity", width = 0.6) +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of Aligned miRNAs per Sample", x = "Case ID", y = "Number of aligned miRNA")
```

```
## Warning: Removed 2 rows containing missing values (geom_bar).
```

![](1_Data_Exploration_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


**Make a new variable bucketing the Processing_time variable into 4 categories: <24 hours, 24-48 hours, 48-72 hours, >72 hours**  
*Hint: use the `case_when` function*  


```r
eDat_plot <- eDat_plot %>% 
  mutate(hours = case_when(
    Processing_time <= 24 ~ "<24 hours",
    between(Processing_time, 24, 48) ~ "24-48 hours", 
    between(Processing_time, 49, 72) ~ "48-72 hours", 
    Processing_time >= 73 ~ ">72 hours",
  ))
```

- Count the number of samples in each of the above formed categories


```r
eDat_plot %>% 
  count (hours) 
```

```
## # A tibble: 4 × 2
##   hours           n
##   <chr>       <int>
## 1 <24 hours   43305
## 2 >72 hours   20209
## 3 24-48 hours 14435
## 4 48-72 hours  8661
```

- A scatter plot showing the expression of hsa-miR-520c-3p per sample for each category, colour by trimester,  `+ facet_wrap(~Trimester)`


```r
eDat_plot %>% 
  mutate(Reads=log2(Reads)) %>% 
  filter(miRNA == "hsa-miR-520c-3p") %>% 
  ggplot( aes(x = hours, y = Reads, color = Trimester )) +
  geom_point() +
  facet_wrap(~Trimester) +
  theme_minimal(base_size = 7) 
```

![](1_Data_Exploration_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


***

### 4.0 Saving files  

**Save `eDat` and `pDat` as `.RDS` objects within your `data` folder


```r
saveRDS(eDat, here::here ("data", "eDat.RDS"))
saveRDS(pDat, here::here ("data", "pDat.RDS"))
```

It was super fun!
***  

