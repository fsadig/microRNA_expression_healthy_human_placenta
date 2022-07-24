---
title: "SSRP: Data Comparison"
author: "Fidan Sadig"
output:
  html_document:
    keep_md: yes
    highlight: haddock #tango, pygments, kate, monochrome, espresso, zenburn, haddock, textmate
    toc: true
    toc_depth: 4
    toc_float:
      collapsed: true
      smooth_scroll: false
    theme: sandstone  #cosmo, paper, lumen, sandstone, simplex, yeti, cerulean, journal, flatly, darkly, readable, spacelab, united
editor_options: 
  chunk_output_type: inline
---  


***  

### **Script #3 Extras** 

***  

### **0.0 Introduction**  

Here, we're going to compare your processed data with the Smith et al., 2020 data, for which you processed it by our method in script 3.  

**Make your own theme and apply to all plots**  


***  

### **1.0 Loading Packages** 


```r
# rmarkdown packages
library(knitr) 
library(rmarkdown)
library(devtools)

# data wrangling packages
library(tidyverse)
library(reshape2)
library(magrittr)
library(stringr)

# data loading packages
library(readxl) 
library(openxlsx)
#library(here)

# plotting/table packages 
library(arsenal)
library(janitor)
library(kableExtra) 
library(egg)
library(RColorBrewer)
library(ggpubr)

# advanced genomic packages
# remotes::install_github("wvictor14/plomics")
library(plomics)
# install.packages("BiocManager")
library(BiocManager)
# BiocManager::install("edgeR")
library(edgeR)
library(pheatmap)

# BiocManager::install(c("limma", "edgeR", "DESeq2", "biomaRt", "miRBaseConverter", "multiMiR", "preprocessCore"))
```


```r
library(sysfonts)
library(showtext)
```

```
## Loading required package: showtextdb
```

```r
library(extrafont)
```

```
## Registering fonts with R
```

```
## 
## Attaching package: 'extrafont'
```

```
## The following object is masked from 'package:showtextdb':
## 
##     font_install
```

```r
my_theme <- theme_minimal() +
  theme(plot.title = element_text(family = "mono", size = 22),
        plot.subtitle = element_text(family = "mono", size = 12),
        legend.text = element_text(family = "mono", size = 10),
        axis.title = element_text(family = "mono", size = 14))


colpal <- list(Trimester = c(`1` = "#cc3399", `2` = "#99FF00", `3` = "#FF3300"),
               Sex = c(`MALE` = "#003366", `FEMALE` = "#FF0066"),
               flowCell = c(`C5JC1ACXX` = "#58508d", `C5JC4ACXX` = "#954693", `C6RGTANXX` = "#cf2d7b"))
```
***


### **2.0 Data**  

Load in your `eDat`, `eNorm_log2`, `pDat`, `mDat`, `smith_eNorm_log2`, `smith_mDat` RDS objects, and your excel sheets containing your DEmiRNAs for each trimester. I added an excel sheet to the data folder called `common_DEmiRNAs.xlsx` which contains the 17 DEmiRNAs that were common between t2_t1, t3_t1, t3_t2  

Also make the `eNorm_plot` object as from before:  



```r
eDat <- readRDS(here::here("data", "eDat.RDS"))
pDat <- readRDS(here::here("data", "pDat.RDS"))
mDat <- readRDS(here::here("data", "mDat.RDS"))
smith_eNorm <-readRDS (here:: here ("data", "smith_eNorm.RDS" ))
smith_mDat <- readRDS(here::here("data", "smith_mDat.RDS"))
eNorm_log2 <-readRDS(here::here("data", "eNorm_log2.RDS"))
  
common_DEmiRNAs <- read_excel(here::here("data", "common_DEmiRNAs.xlsx"))


eNorm_plot <- eNorm_log2 %>% 
  rownames_to_column("mirs") %>% 
  pivot_longer(cols = c(-mirs), names_to = "Case_ID", values_to = "Reads") %>% 
  left_join(pDat, by = "Case_ID") %>% 
  mutate(Trimester = fct_relevel(Trimester, c("1", "2", "3"))) %>% 
  dplyr:: rename("miRNA" = "mirs")
```

### **3.0 EDA on Smith Data**  

First, **go through Melanie's script again and describe her analysis steps** (i.e., what type of filtering did she use, normalization, etc.)   
- She removes samples with the lowest reads 
- simplifies and changes naming of the samples
- converts the table into a dataframe
- removes multimapped miRNA with identical mature sequences
- from metadata, she removes the impure samples
- normalization via TMM (Trimmed Mean of M ) of libraries
- removing biological noise by removing sample with less then 2.5 cmp of their log 2 cpm counts
- pass to removeBatchEffect and remove miRNA with identical sequences



**Is there a correlation between Sex and Oxygenation? Make a table and plot + run and a test**  
*Hint: Which test you would run to compare 2 categorical variables?*   


```r
oxy_sex <- table (smith_mDat$oxygenation, smith_mDat$Sex) %>% 
  chisq.test() 

smith_mDat %>%  
  ggplot (aes(x = oxygenation,  fill = Sex))+
  geom_bar(stat ="count",position = position_dodge2(preserve = "single")) +
  my_theme +
  scale_fill_manual(values = colpal$Sex)+
  labs(title = "Oxygenation vs Sex")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
#X-squared = 2.8393, df = 1, p-value = 0.09198
#results are not significant hence there is not a correlation between Sex and Oxygenation
```


#### **3.1** PCA  

**Do the following:**  

- Run PCA with Trimester, Sex, Oxygenation, and processGroup  
- Make a bar plot of the first 10 PCs
- Make the pval and R-square plots only for the first 10 PCs (*Hint*: Use `dplyr::slice`)
- Plot the PCs by the variables (Make sure the axes are equal)  

**Q1. How many PCs seem to denote the highest variation?**  
PC1 and PC2

**Q2. Which variable seems to most significantly contribute towards the variation in the data?**  
Process group

**Q3. Which 2 variables show an overlap from the PCA plots?**  
Trimester, Sex, Oxygenation


```r
# smith_eLog <- log2(smith_eNorm +1 )

t_eLog <- as.data.frame(t(smith_eNorm))

pca_smith <- prcomp(t_eLog, scale = FALSE, center = TRUE)

scores_sm <- pca_smith$x

summary(pca_smith)
```

```
## Importance of components:
##                            PC1     PC2     PC3     PC4    PC5     PC6     PC7
## Standard deviation     20.8637 14.6961 7.11150 5.70873 5.4719 4.77464 4.34644
## Proportion of Variance  0.3983  0.1976 0.04628 0.02982 0.0274 0.02086 0.01729
## Cumulative Proportion   0.3983  0.5960 0.64227 0.67209 0.6995 0.72035 0.73764
##                            PC8     PC9    PC10    PC11   PC12    PC13    PC14
## Standard deviation     4.01931 3.59637 3.28259 2.98152 2.7849 2.66909 2.61792
## Proportion of Variance 0.01478 0.01184 0.00986 0.00813 0.0071 0.00652 0.00627
## Cumulative Proportion  0.75242 0.76426 0.77412 0.78225 0.7893 0.79587 0.80214
##                           PC15    PC16    PC17    PC18    PC19    PC20    PC21
## Standard deviation     2.53128 2.48385 2.46256 2.36254 2.30451 2.21926 2.20836
## Proportion of Variance 0.00586 0.00565 0.00555 0.00511 0.00486 0.00451 0.00446
## Cumulative Proportion  0.80801 0.81365 0.81920 0.82431 0.82917 0.83368 0.83814
##                           PC22   PC23    PC24   PC25    PC26    PC27    PC28
## Standard deviation     2.17430 2.1178 2.09550 2.0903 2.07901 2.01720 1.97531
## Proportion of Variance 0.00433 0.0041 0.00402 0.0040 0.00396 0.00372 0.00357
## Cumulative Proportion  0.84246 0.8466 0.85059 0.8546 0.85854 0.86226 0.86583
##                           PC29    PC30    PC31    PC32    PC33    PC34    PC35
## Standard deviation     1.95821 1.92978 1.91484 1.88221 1.87575 1.84324 1.81736
## Proportion of Variance 0.00351 0.00341 0.00336 0.00324 0.00322 0.00311 0.00302
## Cumulative Proportion  0.86934 0.87275 0.87611 0.87935 0.88257 0.88568 0.88870
##                           PC36    PC37    PC38    PC39    PC40    PC41   PC42
## Standard deviation     1.80792 1.78839 1.77656 1.75760 1.75647 1.73924 1.7192
## Proportion of Variance 0.00299 0.00293 0.00289 0.00283 0.00282 0.00277 0.0027
## Cumulative Proportion  0.89169 0.89462 0.89751 0.90033 0.90316 0.90593 0.9086
##                           PC43    PC44    PC45    PC46    PC47    PC48    PC49
## Standard deviation     1.71355 1.70040 1.69167 1.66427 1.64299 1.63387 1.60862
## Proportion of Variance 0.00269 0.00265 0.00262 0.00253 0.00247 0.00244 0.00237
## Cumulative Proportion  0.91132 0.91396 0.91658 0.91912 0.92159 0.92403 0.92640
##                           PC50    PC51    PC52    PC53    PC54    PC55    PC56
## Standard deviation     1.60385 1.59265 1.58066 1.57127 1.56688 1.54148 1.52607
## Proportion of Variance 0.00235 0.00232 0.00229 0.00226 0.00225 0.00217 0.00213
## Cumulative Proportion  0.92875 0.93107 0.93336 0.93562 0.93786 0.94004 0.94217
##                           PC57    PC58    PC59   PC60    PC61    PC62   PC63
## Standard deviation     1.51124 1.49246 1.48792 1.4772 1.46792 1.45976 1.4422
## Proportion of Variance 0.00209 0.00204 0.00203 0.0020 0.00197 0.00195 0.0019
## Cumulative Proportion  0.94426 0.94630 0.94832 0.9503 0.95229 0.95424 0.9562
##                           PC64    PC65    PC66    PC67    PC68    PC69    PC70
## Standard deviation     1.43409 1.42126 1.41145 1.39455 1.38410 1.37383 1.36548
## Proportion of Variance 0.00188 0.00185 0.00182 0.00178 0.00175 0.00173 0.00171
## Cumulative Proportion  0.95803 0.95988 0.96170 0.96348 0.96523 0.96696 0.96867
##                           PC71    PC72    PC73    PC74    PC75    PC76    PC77
## Standard deviation     1.34276 1.33619 1.33324 1.32473 1.30790 1.29807 1.28391
## Proportion of Variance 0.00165 0.00163 0.00163 0.00161 0.00157 0.00154 0.00151
## Cumulative Proportion  0.97032 0.97195 0.97358 0.97518 0.97675 0.97829 0.97980
##                           PC78    PC79    PC80    PC81    PC82    PC83    PC84
## Standard deviation     1.27290 1.26620 1.25022 1.24026 1.21444 1.20765 1.20008
## Proportion of Variance 0.00148 0.00147 0.00143 0.00141 0.00135 0.00133 0.00132
## Cumulative Proportion  0.98128 0.98275 0.98418 0.98559 0.98694 0.98827 0.98959
##                           PC85    PC86    PC87    PC88    PC89    PC90    PC91
## Standard deviation     1.18117 1.17451 1.15054 1.13952 1.12802 1.10619 1.09064
## Proportion of Variance 0.00128 0.00126 0.00121 0.00119 0.00116 0.00112 0.00109
## Cumulative Proportion  0.99087 0.99213 0.99334 0.99453 0.99569 0.99681 0.99790
##                           PC92    PC93      PC94
## Standard deviation     1.08184 1.06002 3.172e-14
## Proportion of Variance 0.00107 0.00103 0.000e+00
## Cumulative Proportion  0.99897 1.00000 1.000e+00
```

```r
summary_smith <- data.frame(PC = 1:94, 
                      var_explained = (pca_smith$sdev)^2 / sum((pca_smith$sdev)^2), 
                      cumulative = cumsum(pca_smith$sdev^2 / sum(pca_smith$sdev^2))
                      )

summary_smith <- summary_smith %>% 
  mutate(cumulative_perc = cumulative*100) 
```



```r
cum_bar <- summary_smith %>%
  slice_head(n=10) %>% 
  ggplot(aes(x = PC, y = cumulative_perc)) +
  geom_bar(stat = "identity", fill = "forest green")+
  coord_cartesian(y = c(0,250)) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
  theme_minimal() +
  labs(title = "Cumulative Explained by each PC", x = "Principal Components", y = "Variance Explained", caption = "Normalized Counts")

var_bar<- summary_smith %>%
  slice_head(n=10) %>% 
  ggplot(aes(x = PC, y = var_explained*100)) +
  geom_bar(stat = "identity", fill = "forest green") +
  coord_cartesian(y = c(0,100)) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
  theme_minimal() +
  labs(title = "Variance Explained by each PC", x = "Principal Components", y = "Variance Explained", caption = "Normalized Counts")


ggarrange(cum_bar, var_bar, nrow=1, ncol=2)
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-5-1.png)<!-- -->




```r
smith_mDat_cor <- smith_mDat %>% 
  dplyr::select(Trimester, Sex, oxygenation,processGroup)

smith_mDat_cor <- plomics::lmmatrix(dep = scores_sm, ind = smith_mDat_cor, metric = "Pvalue")


smith_mDat_cor <- t(smith_mDat_cor) 

smith_mDat_cor_10 <- smith_mDat_cor %>% 
  as.data.frame %>% 
  dplyr:: slice_head(n=10)


smith_mDat_cor %>% 
  as.data.frame() %>% 
  mutate(across(everything(), round, 3)) %>% 
  mutate(across(everything(), ~cell_spec(.x, color = ifelse(.x <= 0.05, "green", "")))) %>% 
  kable(escape = F) %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", full_width = F, fixed_thead = T))
```

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Trimester </th>
   <th style="text-align:left;"> Sex </th>
   <th style="text-align:left;"> oxygenation </th>
   <th style="text-align:left;"> processGroup </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.022</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.974</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.027</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.204</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.95</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.482</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.794</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.368</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.009</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.449</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.024</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.447</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.182</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.114</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.975</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.148</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.717</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.875</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.554</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.31</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.464</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.387</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.944</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.278</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.595</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.092</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.744</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.411</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.869</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.73</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.896</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.483</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.681</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.415</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.974</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC11 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.398</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.524</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.875</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.532</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC12 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.682</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.141</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.85</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.695</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC13 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.379</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.157</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.636</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.574</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC14 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.216</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.506</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.858</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.988</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC15 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.992</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.645</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.439</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.881</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC16 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.945</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.888</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.435</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.993</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC17 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.397</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.814</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.105</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.995</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC18 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.752</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.704</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.571</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.805</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC19 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.096</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.859</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.331</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.916</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC20 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.651</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.829</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.519</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.917</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC21 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.914</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.41</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.456</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.96</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC22 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.662</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.502</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.869</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.983</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC23 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.91</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.106</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.849</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.947</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC24 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.575</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.536</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.849</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.762</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC25 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.852</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.203</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.891</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.997</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC26 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.862</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.466</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.629</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.805</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC27 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.548</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.27</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.453</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.9</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC28 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.386</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.128</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.777</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.98</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC29 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.687</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.333</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.614</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.953</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC30 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.7</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.296</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.583</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.918</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC31 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.364</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.09</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.936</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.956</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC32 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.991</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.277</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.923</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.836</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC33 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.665</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.672</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.412</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.935</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC34 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.721</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.012</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.687</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.935</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC35 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.995</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.989</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.969</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.847</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC36 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.538</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.969</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.91</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.931</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC37 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.484</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.626</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.905</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.954</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC38 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.973</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.71</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.754</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.994</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC39 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.985</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.721</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.949</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.825</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC40 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.894</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.084</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.722</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.988</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC41 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.782</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.628</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.743</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.929</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC42 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.727</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.902</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.617</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.969</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC43 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.905</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.02</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.425</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.921</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC44 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.748</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.391</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.894</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.945</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC45 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.554</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.146</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.661</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.952</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC46 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.972</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.059</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.792</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.949</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC47 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.824</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.915</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.548</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.916</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC48 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.798</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.256</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.596</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.947</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC49 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.783</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.812</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.732</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.969</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC50 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.52</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.807</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.62</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.911</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC51 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.803</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.932</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.514</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.941</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC52 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.953</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.663</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.891</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.974</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC53 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.651</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.686</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.761</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.921</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC54 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.947</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.206</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.934</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.889</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC55 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.706</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.708</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.889</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.839</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC56 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.997</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.949</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.46</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.895</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC57 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.793</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.205</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.861</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.99</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC58 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">1</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.012</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.514</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.995</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC59 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.835</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.841</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.531</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.959</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC60 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.537</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.486</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.865</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.986</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC61 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.92</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.143</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.609</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.95</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC62 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.425</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.642</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.87</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.99</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC63 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.443</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.594</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.369</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.932</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC64 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.82</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.92</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.304</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.938</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC65 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.681</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.228</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.575</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.935</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC66 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.704</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.735</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.791</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.981</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC67 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.718</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.925</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.825</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.959</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC68 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.835</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.859</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.977</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.984</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC69 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.719</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.441</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.838</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.94</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC70 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.549</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.742</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.705</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.96</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC71 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.947</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.495</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.787</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.998</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC72 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.797</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.615</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.826</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.988</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC73 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.885</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.372</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.746</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.938</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC74 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.483</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.548</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.887</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.974</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC75 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.675</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.922</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.567</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.86</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC76 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.887</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.134</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.389</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.959</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC77 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.953</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.561</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.801</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.95</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC78 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.782</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.72</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.649</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.935</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC79 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.844</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.212</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.94</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.994</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC80 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.877</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.426</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.589</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.961</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC81 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.419</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.385</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.773</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.876</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC82 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.984</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.809</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.847</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.969</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC83 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.703</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.918</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.939</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.986</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC84 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.739</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.557</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.805</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.998</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC85 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.526</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.423</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.786</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.974</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC86 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.524</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.793</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.828</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.965</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC87 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.692</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.524</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.634</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.937</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC88 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.475</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.39</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.465</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.956</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC89 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.748</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.648</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.846</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.973</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC90 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.76</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.564</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.431</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.984</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC91 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.99</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.328</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.484</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.977</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC92 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.981</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.858</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.385</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.965</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC93 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.414</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.324</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.481</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.972</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC94 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.424</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.011</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.048</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.412</span> </td>
  </tr>
</tbody>
</table>

```r
pc_plot <- smith_mDat_cor_10 %>% 
  mutate(Principal_Component = rownames(smith_mDat_cor_10), PC = as.factor(1:10)) %>% 
  pivot_longer(cols = -c(Principal_Component, PC), names_to = "Variable", values_to = "pval") %>%
  mutate(pval_cat = factor(case_when(
    pval > 0.05  ~ "> 0.05",
    pval < 0.05 & pval > 0.01 ~ "< 0.05",
    pval < 0.01 & pval > 0.001 ~ "< 0.01",
    pval < 0.001 ~ "< 0.001"), 
    levels = c("> 0.05", "< 0.05", "< 0.01", "< 0.001")))

#making a colour palette for the above plot
pc_colpal <- c("white", "#c7dfba", "#8ebe78", "#509e36")
# setting the levels of the pval_cat factor column to the colours in pc_colpal
names(pc_colpal) <- levels(pc_plot$pval_cat)

pval_smith <- pc_plot %>% 
  ggplot(aes(x = PC, y = Variable , fill = pval_cat)) +
  geom_tile(col = "lightgray") +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + #expand function fits the plot to its assigned dimesions 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = pc_colpal)  + 
  coord_fixed() + # very important otherwise height of each tile is elongated!
  labs(title = "Contribution of Variables to Observed Expression\n", x = "", y = "", fill = "P value")
```




```r
pc_rsq_sm <- lmmatrix(dep = scores_sm, ind = smith_mDat_cor, metric = "Rsquared")

pc_rsq_sm_1 <- t(pc_rsq_sm) 

pc_rsq_sm_10 <- pc_rsq_sm_1  %>% 
  as.data.frame %>% 
  slice_head(n=10)

rsq_plot <- pc_rsq_sm_10 %>% 
  mutate(Principal_Component = rownames(pc_rsq_sm_10), PC = as.factor(1:10)) %>% 
  pivot_longer(cols = -c(Principal_Component, PC), names_to = "Variable", values_to = "rsq")

rsq_smith <- rsq_plot %>% 
  ggplot(aes(x = PC, y = Variable , fill = rsq)) +
  geom_tile(col = "lightgray") +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + #expand function fits the plot to its assigned dimesions 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "#509e36") +
  coord_fixed() + # very important otherwise height of each tile is elongated!
  labs(x = "", y = "", fill = "Rsq")


ggarrange(pval_smith, rsq_smith, nrow = 2, ncol = 1) 
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-7-1.png)<!-- -->




```r
scores_sm <- scores_sm %>% 
  as.data.frame() %>% 
  rownames_to_column("samplename")

smith_mDat_pc <- smith_mDat %>% 
  left_join (scores_sm, by = "samplename")

pc1_pc2_1 <- smith_mDat_pc %>% 
  ggplot(aes(x = PC1, y = PC2 , color = Trimester)) +
  geom_point()+
 theme_minimal() +
  coord_cartesian(x = c(-40,40), y = c(-40,40)) +
    labs(title = "Trimester")
  

pc1_pc2_2 <- smith_mDat_pc  %>% 
  ggplot(aes(x = PC1, y = PC2 , color = Sex)) +
  geom_point()+
  theme_minimal() +
   coord_cartesian(x = c(-40,40), y = c(-40,40))+
    labs(title = "Sex")
  
pc1_pc2_3 <- smith_mDat_pc %>% 
  ggplot(aes(x = PC1, y = PC2 , color = oxygenation)) +
  geom_point()+
  theme_minimal() +
  coord_cartesian(x = c(-40,40), y = c(-40,40)) +
    labs(title = "Oxygenation")


pc1_pc2_4<- smith_mDat_pc %>% 
  ggplot(aes(x = PC1, y = PC2 , color = processGroup)) +
  geom_point()+
  theme_minimal() +
  coord_cartesian(x = c(-40,40), y = c(-40,40)) +
  labs(title = "Process Group")

ggarrange(pc1_pc2_1, pc1_pc2_2,pc1_pc2_3, pc1_pc2_4, nrow = 2, ncol = 2)
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

**CHECKPOINT: PC1 should be contributing 40% of the variance**  



### **4.0 Compare DEmiRNA Expression**   

#### **4.1 Robinson DEmiRNAs**  

**Join the 17 DEmiRNAs to their normalized expression**  


```r
eNorm_log2 <- eNorm_log2%>%
  rownames_to_column(var="miRNA") 

common_hm <- common_DEmiRNAs %>% 
  inner_join(eNorm_log2, by = "miRNA")
```

I've added a file called `miRNA_genomic_coordinates.xlsx` which has the details of which chromosome the miRNA lies on, whether it is intronic or exonic, and which is the closest upstream gene to the miRNA if it is intronic.  

**Make a separate df that contains the above information for these 17 miRNAs**  

You'll see that the miRNAs listed are precursors, so **read in the eDat file, and join the precursor names to the mature**



```r
eDat_sp <- eDat %>% 
  rownames_to_column(var="precursor:miRNA") %>% 
  separate("precursor:miRNA", c("precursor", "miRNA"), ":") %>%
  select (precursor:miRNA) %>% 
  distinct

 

common_DEmiRNAs_pr <- common_DEmiRNAs   %>% 
  left_join(eDat_sp ,by = "miRNA" ) 



common_hm_stats <- read_excel(here::here("data", "miRNA_genomic_coordinates.xlsx")) %>% 
  select(mirbase_symbol,chromosome, intragenic_type, distance_from_upstream_exon) %>% 
  dplyr:: rename ("precursor" = "mirbase_symbol" ) %>% 
  inner_join(common_DEmiRNAs_pr, by = "precursor")
```

**Make boxplots of these 17 miRNAs, fill = trimester, and add `facet_grid(cols = vars(chromosome), scale = "free_x")`**


```r
#fig.width=35, fig.height=15
# fig.width=25, fig.height=6
common_hm_stats_eNorm <- common_hm_stats %>% 
  inner_join(eNorm_plot, by = "miRNA")

common_hm_stats_eNorm_a <- common_hm_stats_eNorm %>% 
  filter (chromosome == "chr1" | chromosome == "chr3" | chromosome == "chr5" | chromosome ==  "chr7" | chromosome =="chr9")
  
common_hm_stats_eNorm_b <- common_hm_stats_eNorm %>% 
  filter(chromosome == "chr13" | chromosome == "chr14")

common_hm_stats_eNorm_c<-common_hm_stats_eNorm %>% 
  filter(chromosome == "chr15" |chromosome == "chrX")


gr_a <- common_hm_stats_eNorm_a%>% 
  ggplot(aes (x = miRNA, y= Reads, fill = Trimester))+
  geom_boxplot() +
  my_theme +
  scale_fill_manual(values =  colpal$Trimester)+
  facet_wrap(~miRNA)+
   theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 4))+
  facet_grid(cols = vars(chromosome), scale = "free_x")

gr_b <- common_hm_stats_eNorm_b %>% 
  ggplot(aes (x = miRNA, y= Reads, fill = Trimester))+
  geom_boxplot() +
  my_theme +
  scale_fill_manual(values =  colpal$Trimester)+
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 4))+
  facet_grid(cols = vars(chromosome), scale = "free_x")
  


gr_c <- common_hm_stats_eNorm_c %>% 
  ggplot(aes (x = miRNA, y= Reads, fill = Trimester))+
  geom_boxplot() +
  my_theme +
  scale_fill_manual(values =  colpal$Trimester)+
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 4))+
  facet_grid(cols = vars(chromosome), scale = "free_x",  space = "free")
  


ggarrange(gr_a ,gr_b,gr_c, nrow = 3, ncol = 1, common.legend = TRUE)
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

**Subset to plot the expression of only Chr 14 and Chr X**  


```r
 common_hm_stats_eNorm %>% 
  subset( chromosome == "chr14" | chromosome == "chrX") %>% 
  ggplot(aes (x = miRNA, y= Reads, fill = Trimester))+
  geom_boxplot()+
  my_theme +
  scale_fill_manual(values =  colpal$Trimester)+
  labs(title = "Expression of Common miRNA in Chromosome 14 and Chromosome X") +
  facet_grid(cols = vars(chromosome), scale = "free_x")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


To plot a heatmap to show the differences in the gene expression, we'll first convert our eNorm log values to a [Z-score scale](https://www.simplypsychology.org/z-score.html)


```r
eNorm_log2<- eNorm_log2 %>% 
  column_to_rownames(var="miRNA")

common_hm <- as.data.frame(t(eNorm_log2))
common_hm <- scale(eNorm_log2, center = TRUE, scale = TRUE)
common_hm <- as.data.frame(t(eNorm_log2))
```

**Make a heatmap of the 17 common DEmiRNAs adding Trimester, flowCell, and RQS_cat**  


```r
common_DEmiRNAs_list <- common_DEmiRNAs %>% 
  select (miRNA)

common_DEmiRNAs$miRNA
```

```
##  [1] "hsa-miR-20b-5p"  "hsa-miR-16-5p"   "hsa-miR-106a-5p" "hsa-miR-363-3p" 
##  [5] "hsa-miR-9-3p"    "hsa-miR-590-5p"  "hsa-miR-199b-5p" "hsa-miR-105-5p" 
##  [9] "hsa-miR-431-3p"  "hsa-miR-323a-3p" "hsa-miR-370-3p"  "hsa-miR-432-5p" 
## [13] "hsa-miR-496"     "hsa-miR-767-3p"  "hsa-miR-767-5p"  "hsa-miR-431-5p" 
## [17] "hsa-miR-433-3p"
```

```r
norm_cor_17 <- eNorm_log2 %>% 
  rownames_to_column(var="miRNA") %>% 
  left_join (common_DEmiRNAs_list, by = "miRNA") %>% 
  column_to_rownames(var="miRNA") %>% 
  cor()

pDat <- pDat %>% 
  column_to_rownames(var="Case_ID")

norm_cor_17 %>% 
  pheatmap::pheatmap(clustering_distance_cols = "euclidean", clustering_method = "complete", 
                     cluster_rows = TRUE, show_colnames = TRUE, show_rownames = TRUE,
                     annotation_col = pDat[c("Trimester", "flowCell","RQS")], annotation_row = pDat[c("Trimester", "flowCell","RQS")],
                     main = "17 Common miRNA Sample Correlation")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

Hmm, 4 T2 samples don't cluster together - that's most probably because we have a very low number of genes that we're comparing. If we were to add more genes, all the samples should cluster by trimester (as we saw in our PCA analysis)  


Let's compare the expression of these 17 DEmiRNAs with the Smith data  

**Make a combined plotting df for Robinson and Smith data for norm expression, Trimester, and Sex, and a data column specifying the lab** -  Relevel the data column to order Smith and Robinson  
*Hint: Use thr `rbind` function to bind the Robinson and Smith data*  


```r
smith_eNorm_plot <- smith_eNorm %>% 
  rownames_to_column(var="miRNA") %>% 
  pivot_longer(cols = -c("miRNA"),names_to = "samplename", values_to = "Reads") %>% 
  left_join(smith_mDat, by = "samplename") %>% 
  dplyr:: select (miRNA, Trimester, Sex, Reads) %>% 
   mutate (lab = "Smith") %>% 
   mutate(Trimester = str_replace(Trimester, 'Trimester', "")) %>% 
   mutate(Trimester = fct_relevel(Trimester, c("1", "2")))
```

```
## Warning: Unknown levels in `f`: 1, 2
```

```r
eNorm_plot_cb <- eNorm_plot %>% 
  dplyr:: select(miRNA, Trimester, Sex, Reads) %>% 
   mutate (lab = "Robinson")


combined_plots <-rbind(smith_eNorm_plot, eNorm_plot_cb) %>% 
  mutate (lab = fct_relevel(lab)) 
```

**CHECKPOINT: You should have 239608 rows and 5 columns in your `combined_plots` df**  
I have 238,574 rows and 5 columns

**Filter/Join to only keep the expression for the 17 DEmiRNAs** -  Make sure to use only the trimester 1 and 2 samples, as the Smith study did not include trimester 3 samples.


```r
common_plots <- combined_plots%>% 
  inner_join(common_DEmiRNAs , by = "miRNA") %>% 
  filter (Trimester != 3)
```

Let's try to find if the differences that we see between Robinson and Smith data is significant per trimester  

**Filter to keep either Trimester 1 or 2, x = data, fill = data, use the `stat_compare_means` argument, and move the legend to the bottom of the plot, instead of the default on the right** - stat_compare_means requires specifying the variables (here, our data variable) that you want to compare on the X-axis. 


```r
library(ggpubr)
my_comparisons <- list(c("Robinson","Smith"))
  

 
common_plots %>%
  filter (Trimester != 1) %>% 
  ggplot(aes(x = lab, y = Reads, fill = lab, colour = lab)) +
  geom_boxplot(alpha = 0.8, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#499471", "#ee686d")) +
  scale_colour_manual(values = c("#499471", "#ee686d")) +
  stat_compare_means(comparisons = my_comparisons) +
  coord_cartesian(y = c(0,25)) +
  theme_minimal() +
  facet_wrap(~miRNA) +
  theme(legend.position = "bottom") +
  labs(title = "Comparing Robinson DEmiRNAs with Smith Data: Trimester 2", subtitle = "17 DEmiRNAs")
```

![](4_Data_Comparison_files/figure-html/common plots-1.png)<!-- -->

```r
common_plots %>%
  filter (Trimester != 2) %>% 
  ggplot(aes(x = lab, y = Reads, fill = lab, colour = lab)) +
  geom_boxplot(alpha = 0.8, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#499471", "#ee686d")) +
  scale_colour_manual(values = c("#499471", "#ee686d")) +
  stat_compare_means(comparisons = my_comparisons) +
  coord_cartesian(y = c(0,25)) +
  theme_minimal() +
  facet_wrap(~miRNA) +
  theme(legend.position = "bottom") +
  labs(title = "Comparing Robinson DEmiRNAs with Smith Data: Trimester 1", subtitle = "17 DEmiRNAs")
```

![](4_Data_Comparison_files/figure-html/common plots-2.png)<!-- -->

```r
 common_plots %>%
  filter (Trimester != 3) %>% 
  ggplot(aes(x = lab, y = Reads, fill = lab, colour = lab)) +
  geom_boxplot(alpha = 0.8, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#499471", "#ee686d")) +
  scale_colour_manual(values = c("#499471", "#ee686d")) +
  stat_compare_means(comparisons = my_comparisons) +
  coord_cartesian(y = c(0,25)) +
  theme_minimal() +
  facet_wrap(~miRNA, ncol = 3) +
  theme(legend.position = "bottom") +
  labs(title = "Comparing Robinson DEmiRNAs with Smith Data: Trimester 1 and 2 ", subtitle = "17 DEmiRNAs")
```

![](4_Data_Comparison_files/figure-html/common plots-3.png)<!-- -->

**Based on the above two plots and the pvalues, which miRNAs show similar expression across both studies?**  
hsa-miR-370-3p with a pvalue of 0.29 in trimester 1 and 2.

In trimester 2, hsa-miR-363-3p,hsa-miR-433-3p, hsa-miR-767-3p.

In trimester 1, hsa-miR-105-5p,hsa-miR-767-5p,hsa-miR-370-3p.

#### **4.2 Smith DEmiRNAs** 

**Compare the normalized expression of the top 10 DEmiRNAs found in the Smith study using `combined_plots`** - Make the same plots as above, separate for both Trimester 1 and 2. 


```r
smith_de <- c(
"hsa-miR-30d-5p",
"hsa-miR-125a-5p",
"hsa-miR-517a-3p",
"hsa-miR-199a-3p",
"hsa-miR-26b-5p",
"hsa-miR-26a-5p",
"hsa-let-7a-5p",
"hsa-miR-21-5p",
"hsa-miR-126-3p",
"hsa-miR-516b-5p"
)


combined_plots %>% 
  filter (miRNA == smith_de) %>% 
  filter (Trimester != 1) %>% 
  ggplot(aes(x = lab, y = Reads, fill = lab, colour = lab)) +
  geom_boxplot(alpha = 0.8, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#499471", "#ee686d")) +
  scale_colour_manual(values = c("#499471", "#ee686d")) +
  stat_compare_means(comparisons = my_comparisons) +
  coord_cartesian(y = c(0,25)) +
  theme_minimal() +
  facet_wrap(~miRNA) +
  theme(legend.position = "bottom") +
  labs(title = "Comparing top 10 Smith DEmiRNA: Trimester 2", subtitle = "17 DEmiRNAs")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
 combined_plots %>% 
  filter (miRNA == smith_de) %>% 
  filter (Trimester != 2) %>% 
  ggplot(aes(x = lab, y = Reads, fill = lab, colour = lab)) +
  geom_boxplot(alpha = 0.8, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#499471", "#ee686d")) +
  scale_colour_manual(values = c("#499471", "#ee686d")) +
  stat_compare_means(comparisons = my_comparisons) +
  coord_cartesian(y = c(0,25)) +
  theme_minimal() +
  facet_wrap(~miRNA) +
  theme(legend.position = "bottom") +
  labs(title = "Comparing top 10 Smith DEmiRNA: Trimester 1", subtitle = "17 DEmiRNAs")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

```r
 combined_plots %>% 
  filter (miRNA == smith_de) %>% 
  filter (Trimester != 3) %>% 
  ggplot(aes(x = lab, y = Reads, fill = lab, colour = lab)) +
  geom_boxplot(alpha = 0.8, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#499471", "#ee686d")) +
  scale_colour_manual(values = c("#499471", "#ee686d")) +
  stat_compare_means(comparisons = my_comparisons) +
  coord_cartesian(y = c(0,25)) +
  theme_minimal() +
  facet_wrap(~miRNA) +
  theme(legend.position = "bottom") +
  labs(title = "Comparing top 10 Smith DEmiRNA: Trimester 1 and 2", subtitle = "17 DEmiRNAs")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-17-3.png)<!-- -->


**Based on the above two plots, which miRNAs show similar expression across both studies?**  
Trimester  2  hsa-miR-7a-5p 

**Select only the one miRNA which shows no diff in exp and redo the plot** - add the argument `caption` within `labs` and write down how many samples are present for each study


```r
common_plots %>%
  filter (Trimester != 2) %>%
  filter (miRNA  == "hsa-miR-370-3p") %>% 
  ggplot(aes(x = lab, y = Reads, fill = lab, colour = lab)) +
  geom_boxplot(alpha = 0.8, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#499471", "#ee686d")) +
  scale_colour_manual(values = c("#499471", "#ee686d")) +
  stat_compare_means(comparisons = my_comparisons) +
  coord_cartesian(y = c(0,25)) +
  theme_minimal() +
  facet_wrap(~miRNA) +
  theme(legend.position = "bottom") +
  labs(title = " DEmiRNAs with no change in expression", caption = "Smith: 94, Robinson:40")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
common_plots %>% 
  filter(Trimester == 1) %>% 
  ggplot(aes(x = lab, y = Reads, fill = lab, colour = lab)) +
  geom_boxplot(alpha = 0.8, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#499471", "#ee686d")) +
  scale_colour_manual(values = c("#499471", "#ee686d")) +
  stat_compare_means(comparisons = my_comparisons) + 
  coord_cartesian(y = c(0,25)) +
  theme_minimal() +
  facet_wrap(~miRNA) +
  theme(legend.position = "bottom") +
  labs(title = "Comparing Robinson DEmiRNAs with Smith Data: Trimester 1", subtitle = "17 DEmiRNAs")
```

```
## Warning: Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
## Computation failed in `stat_signif()`:
## not enough 'y' observations
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-18-2.png)<!-- -->

#### **4.3 miRNA Clusters**

I've added an excel file that has the 3 major miRNA clusters showing expression in the placenta - we'll just be looking at C14MC and C19MC  

**Separate the clusters to make new dfs of the respective cluster** 


```r
clusters <- read_excel(here::here("data", "miRNA_Clusters.xlsx"))

c14 <- clusters %>% 
 filter(grepl ("14", cluster) )

c19 <-  clusters %>% 
 filter(grepl ("19", cluster) )
```

We're going to check how much of the total miRNA expression we see is by the clusters, i.e., the proportion of expression of these clusters within each of our samples 

**Get the total reads for all miRNAs per sample within eNorm**  
*Hint* - Use the `colsums` function


```r
colsums_roblab <- eNorm_log2 %>% 
  colSums (.) %>% 
  as.data.frame()
```

**Get the total reads for the C14MC and C19MC miRNAs per sample**


```r
eNorm_log2 <- eNorm_log2 %>% 
  rownames_to_column("miRNA")

c14_mirs <- c14 %>% 
  dplyr::select(miRNA) %>% 
  left_join(eNorm_log2, by = "miRNA") 
```

You'll see that the df is mostly empty - can you guess why?  

We'll make a separate df from eNorm which removes the `-3p` and `-5p` from the end of the miRNA - making them the precursor version of the miRNA.  

We'll use the `str_replace_all` function:  


```r
eNorm_pre <- eNorm_log2 

eNorm_pre$miRNA <-  str_replace_all(eNorm_pre$miRNA , pattern = "-5p", replacement = "") 
eNorm_pre$miRNA <-  str_replace_all(eNorm_pre$miRNA , pattern = "-3p", replacement = "") 

c14_mirs <- c14 %>% 
  dplyr::select(miRNA) %>% 
  left_join(eNorm_pre, by = "miRNA") %>% 
  rename("precursor" = "miRNA")

# removing the rows with NA
c14_mirs <- na.omit(c14_mirs) 


#do the same for c19
c19_mirs <- c19 %>% 
  dplyr::select(miRNA) %>% 
  left_join(eNorm_pre, by = "miRNA")%>% 
  rename("precursor" = "miRNA")


c19_mirs <- na.omit(c19_mirs)
```

**Get the total reads for both c14 and c19**  


```r
# c14_colsums <- as.data.frame(colSums(c14_mirs))
# 
# c19_colsums <- as.data.frame(colSums(c19_mirs))
```

You should get an error here: `Error in colSums(c14_prop) : 'x' must be numeric`. That's because we have our miRNAs in a column, so our df is not completely a numerical one. 

**Go ahead and try converting miRNAs into rownames**  


```r
# c19_mirs <- c19_mirs %>%  
#   as.data.frame() %>% 
#   column_to_rownames("miRNA")
```


You'll get an error stating `Error in .rowNamesDF<-(x, value = value) :  duplicate 'row.names' are not allowed` -- and that's because we converted our mature miRNAs to precursor to be able to join the cluster and eNorm dfds above.  

**Get the colsums of c14 and c19 by using data manipulation functions of your choosing**  


```r
list_mirna_14 <- c14_mirs$precursor

colsums_c14 <- c14_mirs %>% 
  select (-precursor) %>% 
  colSums(.) %>% 
  as.data.frame()


list_mirna_19 <- c19_mirs$precursor

colsums_c19 <- c19_mirs %>% 
  select (-precursor) %>% 
  colSums(.) %>% 
  as.data.frame()
```

**Create a  stacked bar plot for each sample showing the total number of reads, and the proportion of C14 and C19 within them by percentage** - The total bar chart should add up to a 100 (%)


```r
all_colsums <- cbind(colsums_roblab, colsums_c14, colsums_c19) 
colnames(all_colsums) <- c("Total", "C14MC", "C19MC") 

all_colsums_C19 <- all_colsums %>% 
  rownames_to_column(var="samples")  %>% 
  mutate(Clusters = "C19_p")


all_colsums_C14 <- all_colsums %>% 
  rownames_to_column(var="samples")  %>% 
  mutate(Clusters = "C14_p")

all_colsums_total <- all_colsums %>% 
   rownames_to_column(var="samples")  %>% 
   mutate(Clusters = "rest_p")


all_colsums_up <- bind_rows (all_colsums_C19 ,all_colsums_C14,all_colsums_total)


all_colsums_up <- all_colsums_up %>% 
  mutate(Percentage = case_when(
    Clusters == "C19_p" ~ C19MC*100/Total,
     Clusters =="C14_p" ~ C14MC*100/Total,
     Clusters == "rest_p" ~ (Total-C14MC-C19MC)*100/Total))

   
all_colsums_up %>%
  ggplot( aes( x = samples, y = Percentage, fill = Clusters )) +
    geom_bar(position="fill", stat="identity")+
    scale_y_continuous(labels = scales::percent_format())
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

**Let's plot all the miRNAs within each cluster by trimester and also by sex**

However, now that we have the precursor, all our mature miRNAs will be grouped into one miRNA. And we can't just add a -3p and -5p to the ends of the precursors, because some miRNAs don't have a 3p or 5p version (they just have the one version)  

**Join the precursors to the mature miRNA names**  

You'll see that when you try to join, you don't get any matches, but when you view the df and search, the miRNAs match - can you figure out why the df returned has no matches?  

**Fix the error on why the dfs aren't joining**  
*Hint: Look at the spelling of the miRNAs carefully. And it's a new function that I've introduced in this script somewhere above*    



```r
eDat_up <-  eDat %>% 
  rownames_to_column(var = "pre:miRNA") 

pre_mature <-eDat_up %>% 
  select ("pre:miRNA") %>% 
   separate("pre:miRNA", c("precursor", "miRNA"), ":") 

pre_mature$precursor <-  str_replace_all(pre_mature$precursor   , pattern = "mir", replacement = "miR") 

#C19 cluster
c19_mirs$precursor  <- str_replace_all(c19_mirs$precursor   , pattern = "mir", replacement = "miR") 

c19_mirs_joined <- c19_mirs %>% 
    left_join(pre_mature, by = "precursor")  %>% 
    pivot_longer(cols = -c("miRNA", "precursor"),names_to = "Case_ID", values_to = "Reads")

#C14 cluster

c14_mirs$precursor  <- str_replace_all(c14_mirs$precursor   , pattern = "mir", replacement = "miR") 

c14_mirs_joined <- c14_mirs %>% 
    left_join(pre_mature, by = "precursor")  %>% 
    pivot_longer(cols = -c("miRNA", "precursor"),names_to = "Case_ID", values_to = "Reads")
```

**Plot the C14MC and C19MC clusters by trimester/sex and facet by miRNA**  


```r
trimester_sex <- eNorm_plot %>% 
  select(miRNA, Trimester, Sex)

c19_mirs_joined_eNorm <- c19_mirs_joined %>% 
  left_join (trimester_sex, by = "miRNA")

c19_mirs_joined<- c19_mirs_joined %>% 
  left_join (mDat, by = "Case_ID")  #when combined with the timester_sex  df I got NA

c14_mirs_joined <- c14_mirs_joined %>% 
  left_join (mDat, by = "Case_ID")

c19_mirs_joined %>% 
  dplyr::count(Trimester = NA)
```

```
## # A tibble: 1 × 2
##   Trimester     n
##   <lgl>     <int>
## 1 NA         2910
```

```r
c19_mirs_joined%>% 
  # filter (!is.na(c19_mirs_joined$miRNA)) %>% #filtering NA that I got from trimester_sex df
  ggplot (aes( y = Reads, color = Trimester))+
  geom_boxplot() +
  my_theme +
  scale_fill_manual(colpal$Trimester) +
  coord_cartesian(y = c(0,20))+
  facet_wrap(~miRNA) +
  labs(title = " Expression of microRNAs in C19 cluster by Trimester")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

```r
c19_mirs_joined%>% 
  ggplot (aes( y = Reads, color = Sex))+
  geom_boxplot() +
  my_theme +
  scale_fill_manual(colpal$Sex) +
  coord_cartesian(y = c(0,20))+
  facet_wrap(~miRNA) +
  labs(title = " Expression of microRNAs in C19 cluster by Sex")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-28-2.png)<!-- -->

```r
c14_mirs_joined%>% 
  ggplot (aes( y = Reads, color = Trimester))+
  geom_boxplot() +
  my_theme +
  scale_fill_manual(colpal$Trimester) +
  coord_cartesian(y = c(0,20))+
  facet_wrap(~miRNA) +
  labs(title = " Expression of microRNAs in C14 cluster by Trimester")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-28-3.png)<!-- -->

```r
c14_mirs_joined%>% 
  ggplot (aes( y = Reads, color = Sex))+
  geom_boxplot() +
  my_theme +
  scale_fill_manual(colpal$Sex) +
  coord_cartesian(y = c(0,20))+
  facet_wrap(~miRNA) +
  labs(title = " Expression of microRNAs in C14 cluster by Sex")
```

![](4_Data_Comparison_files/figure-html/unnamed-chunk-28-4.png)<!-- -->

***  

> You're done with your project!  

***  
