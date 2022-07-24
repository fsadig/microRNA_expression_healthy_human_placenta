---
title: "SSRP: Quality Control"
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

---


***  

### **Script #1 Extras**   

Use this space to add additionally notes or plots for yourself from what we discussed after going through script #1  


```r
# data loading packages
library(readxl) 
library(openxlsx)
library(here)
```

```
## here() starts at /Users/fidan/BCCHR/SSRP
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
## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
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
#I had to put these lines here too to run the plots 

eDat <- readRDS(here::here("data", "eDat.RDS"))  
pDat <- readRDS(here::here("data", "pDat.RDS"))


eDat_nonames <- eDat %>% 
  rownames_to_column(var="Precursor:miRNA") %>% 
  separate("Precursor:miRNA", c("precursor", "miRNA"), ":")


eDat_plot <- eDat_nonames %>% 
  pivot_longer(cols = -c(precursor, miRNA), names_to = "Case_ID", values_to = "Reads")


eDat_plot <- eDat_plot %>% 
  inner_join( pDat, c ("Case_ID"))


eDat_plot %>% 
  mutate(Reads=log2(Reads)) %>% 
  filter (miRNA == c("hsa-miR-21-5p","hsa-miR-200a-5p","hsa-miR-373-5p")) %>% 
  ggplot(aes(x = miRNA, y = Reads, fill = Trimester)) +
  geom_boxplot() +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of Aligned miRNAs per Sample", x = "Case ID", y = "Number of aligned miRNA")
```

```
## Warning: Removed 2 rows containing non-finite values (stat_boxplot).
```

![](2_QC_files/figure-html/boxplot-1.png)<!-- -->




Code chunk parameters:  
- eval: whether you want to evaluate the output of that chunk (sometimes you'll write a few different code versions performing the same
functions. In order to keep them separate, I put them into different chunks and set the ones I don't want to use to eval = FALSE)  
- echo: whether you want to output the code in that chunk (I set my package loading chunk to echo = FALSE)  
- warning: display any warnings associated with the code ran in that chunk  
- message: display any messages associated with the code ran in that chunk  
- error: display any errors associated with the code ran in that chunk  
- results: options to display the code output  

Check the following on how to set custom themes for your plots:  



```r
#removed eval = False to run the my_theme

# load any font from Google Fonts that you'd like - my personal favourite is Montserrat. Check out all the fonts here: https://fonts.google.com

# install.packages("sysfonts")
# install.packages("showtext")
# install.packages("extrafont")

library(sysfonts)
library(showtext)
library(extrafont)

font_add_google("Montserrat", "Mont") # I'm saying to add the Montserrat font, and name it as Mont

showtext_auto() # required to render custom fonts
knitr::opts_chunk$set(fig.showtext = TRUE, fig_retina = 1)

# my custom ggplot theme
my_theme <- theme_minimal() +
  theme(plot.title = element_text(family = "Mont", size = 14),
        plot.subtitle = element_text(family = "Mont", size = 12),
        legend.text = element_text(family = "Mont", size = 10),
        axis.title = element_text(family = "Mont", size = 14))

# my custom colour palette for my variables - I like to use this website: https://www.schemecolor.com/ 
# here, we're making another data type called a list - a list can stire multiple dataframes. We're making a separate df for each our pur variables and assigning it a corresponding hex code (colour code) which will then be output when we specify the manual colour/fill option in ggplot

colpal <- list(Trimester = c(`1` = "#ecb3da", `2` = "#c36cac", `3` = "#981580"),
               Sex = c(`MALE` = "#E05D5D", `FEMALE` = "#ffb708"),
               flowCell = c(`C5JC1ACXX` = "#58508d", `C5JC4ACXX` = "#954693", `C6RGTANXX` = "#cf2d7b"))

# add however many more variables that you'd like!

# when you want to specify these pariculat colours, within your ggplot add: 
 
scale_colour_manual(values = colpal$whichever_variable_you_want_to_colour) 
# OR 
scale_fill_manual(values = colpal$whichever_variable_you_want_to_fill)
```


### **0.0 Introduction**   

Here, we're going to perform quality control steps on our data.  

What is QC and why is it important?  

Before we start analysing our data and getting meaningful interpretations of our result, we have to ensure that the data that we're actually using is representative of our sample expression. We don't want samples or sequences (here, miRNAs) to represent a biased representation. For example, the majority of sequences/species will record 0 counts (transcripts), and keeping these sequences in our data would lead us to believe that the average expression of our whole dataset is close to zero (when that is not the case). And hence, the basic steps of QC include:  

- Filtering: removing sequences/samples of low quality  
- Normalization: Re-adjusting and standardising expression of filtered sequences for each gene between samples, and further all genes within the same sample  
- Hierarchical Clustering: Identifying the variation between samples (how similar or disimilar are they from each other)  
- Dimensionality Reduction: Identifying the variables within our data that contribute to the variation observed in our data   

[This workflow](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html) provides good explanations of the steps that we're going to be using.  


***  

### **1.0 Loading Packages and Data** 


```r
# rmarkdown packages
library(knitr) 
library(rmarkdown)
library(devtools)



# plotting/table packages 
library(arsenal)
library(janitor)
library(kableExtra) 
library(egg)
library(RColorBrewer)
library(ggpubr)

# advanced genomic packages
remotes::install_github("wvictor14/plomics")
library(plomics)
#install.packages("BiocManager")
library(BiocManager)
BiocManager::install("edgeR")
library(edgeR)
library(pheatmap)

# BiocManager::install(c("limma", "edgeR", "DESeq2", "biomaRt", "miRBaseConverter", "multiMiR", "preprocessCore"))
```

***

### **2.0 Data**  

Load in your `eDat` and `pDat` RDS (R objects) saved from the `Data_Exploration` script 


```r
eDat <- readRDS(here::here("data", "eDat.RDS"))
pDat <- readRDS(here::here("data", "pDat.RDS"))
```



Okay, first we're going to visualise the spread of our data. For this, we'll use something called as a density plot, as well as boxplots.  

Go back to the `Data_Exploration` script, and find the `pivot_longer` function. We'll be using that to convert eDat into a longer df to be able to plot.  

**Convert eDat into a longer df, with `names_to` to Case_ID and `values_to` to Reads, and join pDat by Case_ID**  

Convert the rownames to a column first


```r
eDat_plot <- eDat %>% 
  rownames_to_column(var = "miRNA") %>% 
  pivot_longer(cols = -c("miRNA"),names_to = "Case_ID", values_to = "Reads") %>% 
  inner_join(., pDat, by = "Case_ID")
```


```r
eDat_plot %>% 
  ggplot(aes(x = Reads, colour = Case_ID)) +
  geom_density() +
  my_theme
```

![](2_QC_files/figure-html/raw-spread-1.png)<!-- -->

Hmm, so that doesn't show us much - and the reason being that the spread of our data is too wide, and hence we're not able to visualise it well. As most of the expression it at the extremes, with the majority of sequences showing wither 0 or very high expression, we're not able to see the spread.  

To overcome this, we're going to convert our expression values to a log10 scale (expression data is almost always visualised by converting to a log scale (usually log2 or log10) to be able to visualise the relative differences)


```r
eDat_plot %>% 
  ggplot(aes(x = log10(Reads), colour = Case_ID)) +
  geom_density() +
  my_theme +
  theme(legend.position = "none") #to remove the legends of our plot
```

```
## Warning: Removed 41558 rows containing non-finite values (stat_density).
```

![](2_QC_files/figure-html/raw-spread-log-1.png)<!-- -->


However, we get a warning saying `Warning: Removed 41558 rows containing non-finite values (stat_density).` Can you guess why?  

Go ahead and convert 0 to its log10 expression


```r
log10(0)
```

```
## [1] -Inf
```

And here we see that the [log2 value of 0 is infinity](https://www.google.com/search?q=log2+value+of+0&rlz=1C1CHBD_enCA1009CA1009&oq=log2+value+of+0&aqs=chrome..69i57j0i512j0i22i30j0i390l4.382j0j7&sourceid=chrome&ie=UTF-8), and we cannot plot an infinite value.  

And so, to be able to plot 0 values, we do a `log2(x+1)` transformation, the most common way to be able to visualise these sequences. This transformation has a pseudo count/value of 1 to all values, essentially shifting our entire expression matrix by +1. 


```r
eDat_plot %>% 
  ggplot(aes(x = log10(Reads + 1), colour = Case_ID)) +
  geom_density() +
  theme_minimal() +
  theme(legend.position = "none")
```

![](2_QC_files/figure-html/raw-spread-log2-1-1.png)<!-- -->

And this shows us the spread of our raw data. What do you think is the interpretation of this plot?  

It seems like it is left skewed. It shows that majority of reads are in the 0-1 fold range. 

**Facet the above density plot by trimester**  


```r
eDat_plot %>% 
  ggplot(aes(x = log10(Reads + 1), colour = Case_ID)) +
  geom_density() +
  theme_minimal() +
  facet_wrap("Trimester")+
  theme(legend.position = "none")
```

![](2_QC_files/figure-html/density-trimester-1.png)<!-- -->


Let's now see the expression spread per sample using boxplots:


```r
eDat_plot %>% 
  ggplot(aes(x = Case_ID, y = log10(Reads + 1), fill = Case_ID)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") 
```

![](2_QC_files/figure-html/raw-boxplots-1.png)<!-- -->


We're now going to assign the plots to a variable, and display them together:  

**Assign both plots to variables called d1 (density) and b1 (boxplot 1) and output them together**  

Add `+ labs(title = "Gene Expression", subtitle = "Raw Counts")` to both plots, and `+ theme(axis.text.x = element_text(angle = 35, hjust = 1))` to the boxplot  


```r
d1 <- eDat_plot %>% 
  ggplot(aes(x = log10(Reads + 1), colour = Case_ID)) +
  geom_density() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Gene Expression", subtitle = "Raw Counts")

b1 <-  eDat_plot %>% 
  ggplot(aes(x = Case_ID, y = log10(Reads + 1), fill = Case_ID)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none")  +
  labs(title = "Gene Expression", subtitle = "Raw Counts") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))


ggarrange(d1, b1)
```

![](2_QC_files/figure-html/raw-plots-1.png)<!-- -->

***

### **3.0 Filtering**  

Usually, filtering is performed by keeping a certain number of reads per miRNA in a certain number of samples. So, for example:  

- By Counts: Expression Count of 1 or more in all samples  
- By Counts and % of Samples: Expression Count of 1 or more in 10% of samples  
- By Counts and No. of Samples: Expression Count of 2 or more in 2 samples  

#### **3.1 Removing miRNAs with 0 counts**  




```r
eFil_0 <- eDat %>% 
  filter_all(any_vars(. > 0))
```

**To check if we have any miRNAs with 0 in all samples mistakenly left, use the `RowSums` function to make a new column and `filter` to select entries which have an entry of 0 within the sum column**  


```r
eFil_0 %>% 
  mutate(Sums = rowSums(.)) %>% 
  filter (Sums == 0) 
```

```
##  [1] FT42  FT52  FT54  FT60  FT67  FT73  FT85  FT65  FT82  FT74  PL166 PL165
## [13] PL167 FT36  FT41  PM92  PL149 PL170 FT91  PM84  PL148 PL171 PM114 PM120
## [25] PM167 PM181 PM135 PM153 PM112 NTD4  Sums 
## <0 rows> (or 0-length row.names)
```

```r
 eFil_0 %>% 
  mutate(Sums = rowSums(.)) %>% 
  glimpse() 
```

```
## Rows: 2,548
## Columns: 31
## $ FT42  <dbl> 444, 96079, 3, 96112, 450, 96864, 43, 7718, 3, 6968, 226, 137, 1…
## $ FT52  <dbl> 647, 164316, 2, 164372, 649, 165176, 78, 12420, 35, 14281, 476, …
## $ FT54  <dbl> 709, 293191, 22, 293225, 712, 294468, 294, 46712, 19, 39354, 828…
## $ FT60  <dbl> 813, 226991, 12, 227042, 823, 228053, 127, 31511, 74, 12717, 981…
## $ FT67  <dbl> 495, 186434, 8, 186507, 501, 187442, 66, 12363, 18, 23337, 292, …
## $ FT73  <dbl> 314, 55365, 3, 55392, 314, 55763, 25, 4229, 8, 6943, 161, 134, 1…
## $ FT85  <dbl> 276, 63794, 0, 63813, 278, 64435, 37, 7375, 7, 3376, 161, 99, 11…
## $ FT65  <dbl> 1493, 522367, 16, 522475, 1506, 524275, 214, 49933, 38, 20605, 1…
## $ FT82  <dbl> 1183, 396173, 3, 396201, 1196, 398205, 80, 25687, 33, 17536, 583…
## $ FT74  <dbl> 1483, 306320, 7, 306488, 1498, 308350, 160, 27529, 50, 13747, 98…
## $ PL166 <dbl> 37, 4892, 2, 4885, 37, 4904, 1, 573, 0, 763, 24, 49, 15, 741, 1,…
## $ PL165 <dbl> 377, 34289, 2, 34249, 380, 34388, 10, 4663, 13, 5342, 199, 272, …
## $ PL167 <dbl> 26, 4281, 0, 4278, 26, 4291, 3, 1333, 3, 972, 37, 65, 15, 1172, …
## $ FT36  <dbl> 82, 18877, 1, 18851, 82, 18994, 8, 2981, 3, 2978, 73, 169, 18, 3…
## $ FT41  <dbl> 161, 54184, 1, 54159, 163, 54412, 10, 5212, 0, 7149, 81, 226, 40…
## $ PM92  <dbl> 205, 53078, 5, 53113, 205, 53148, 49, 6631, 6, 3194, 167, 191, 6…
## $ PL149 <dbl> 373, 97396, 2, 97363, 374, 97813, 40, 6804, 3, 6575, 267, 364, 1…
## $ PL170 <dbl> 63, 8462, 1, 8458, 63, 8505, 1, 643, 5, 568, 25, 40, 27, 1884, 0…
## $ FT91  <dbl> 284, 36696, 0, 36658, 287, 36976, 21, 2567, 4, 2243, 190, 175, 9…
## $ PM84  <dbl> 801, 171395, 3, 171415, 803, 171989, 176, 39280, 18, 11563, 1017…
## $ PL148 <dbl> 717, 154473, 9, 154346, 723, 155008, 51, 8912, 18, 9547, 351, 50…
## $ PL171 <dbl> 43, 3864, 3, 3867, 44, 3871, 2, 524, 3, 262, 22, 24, 20, 1374, 1…
## $ PM114 <dbl> 1761, 382620, 5, 382570, 1765, 383478, 202, 39757, 48, 24359, 16…
## $ PM120 <dbl> 146, 31237, 0, 31180, 146, 31221, 20, 2763, 6, 1372, 70, 144, 55…
## $ PM167 <dbl> 737, 186118, 8, 185711, 745, 185863, 74, 13911, 15, 8704, 412, 9…
## $ PM181 <dbl> 120, 27350, 0, 27281, 122, 27375, 7, 2578, 2, 1433, 77, 158, 50,…
## $ PM135 <dbl> 845, 161419, 5, 161003, 851, 161174, 54, 14197, 22, 9396, 557, 6…
## $ PM153 <dbl> 613, 109686, 13, 109434, 614, 109331, 70, 13280, 21, 6478, 568, …
## $ PM112 <dbl> 1103, 214843, 4, 214606, 1104, 215051, 201, 29726, 26, 11124, 95…
## $ NTD4  <dbl> 944, 186632, 4, 186168, 959, 187011, 110, 12151, 29, 10555, 490,…
## $ Sums  <dbl> 17295, 4252822, 144, 4251222, 17420, 4267834, 2234, 433963, 530,…
```

**CHECKPOINT:** You `eFil_0` df should have 2548 rows and 31 columns  

Awesome!

However, that does leave us with some miRNAs that have expression > 0 in some samples, but with expression = 0 in the rest of the samples  

**Separate column 1 in eDat_plot, and pull out hsa-miR-10226 and make a bar chart of it's expression by Sample**  


```r
eDat_plot %>% 
  separate("miRNA", c("precursor", "miRNA"), ":") %>% 
  filter (miRNA == "hsa-miR-10226") %>% 
  ggplot(aes(x = Case_ID , y = Reads )) +
  geom_bar(stat = "identity", width = 0.6) +
  my_theme +
  labs(title = "Expression of hsa-miR-10226 by sample", x = "Sample", y = "Reads")
```

![](2_QC_files/figure-html/plot_0-1.png)<!-- -->

Do you think we should filter out miRNAs with 0 expression in all samples, or only miRNAs with 0 expression in some samples? Give a pro and con that you can think for both filtering criteria:  

Filter out miRNAs with 0 expression in all samples:  
Pro: higher average per sample and a better representation of the landscape. 
Con: loose information about microRNA that do not contribute at all which can be useful for future researches.

Filter out miRNAs with 0 expression in some samples:  
Pro: decreasing noise and focusing on the significant samples. 
Con:loosing potential candidates to include in the final result. microRNAs maybe have the zero in some samples because wrong extractions or some other false signals???

#### **3.2 Filtering by % of Samples**     

Usually filtering of % of samples occurs when you have distict groups, i.e., cases and controls, for examples. You should select to keep miRNAs with a count of x in x% of samples in each of the groups, and then retain the miRNAs that pass that filtering criteria. Some of the miRNAs that pass would be common between the two groups, and some would be present in only either of the two groups  
Here, filtering by % of samples doesn't make the most sense, as even if we say 10% of samples, that 10 samples (since we have a total of 30 samples) and those 10 samples 


#### **3.3 Filtering by Number of Samples**      

Let's check the number of samples we have per group by 2 of the variables that we can hypothesize will be contributing towards differences in expression - trimester and sex  


```r
pDat %>% 
  tabyl(Trimester, Sex)
```

```
##  Trimester FEMALE MALE
##          1      3    2
##          2     10    6
##          3      4    5
```

Okay, so the group with the lowest number of samples is Trimester 1 Males - 2 samples. If we were to think about the genes that might be differentially expressed, we can assume that there might be a few that show differential expression in these two samples as compared to the rest. Now, 2 does seem to be too low of a number to filter by, but given that we have only 30 samples, it's fine.

So, we're going to filter to only keep miRNAs with an expression count of 1 or more in 2 or more samples 

For that, we'll first convert our `eFil_0` into a logical df, giving us a TRUE value is expression is equal to or greater than 1 count, FALSE is expression is 0, and then sum the number of TRUEs we have, and if the sum is more than 2, then those are the miRNAs we keep


```r
eFil_1_2 <- eFil_0 >=1

eFil_1_2 <- eFil_1_2 %>% 
  as.data.frame() %>% #the above function converts the df into a matrix, and there are a few functions we cannot apply to matrices
  mutate(Sum = rowSums(.))

eFil_1_2 <- eFil_1_2 %>% 
  filter(Sum >= 2)

eFil_0 <- eFil_0 %>% 
  rownames_to_column("mirs")
  
eFil <- eFil_1_2 %>% 
  rownames_to_column("mirs") %>% 
  dplyr::select(mirs) %>% 
  left_join(eFil_0, by = "mirs") 

all(eFil$mirna == rownames(eFil_1_2))
```

```
## [1] TRUE
```

```r
glimpse(eFil)
```

```
## Rows: 2,344
## Columns: 31
## $ mirs  <chr> "hsa-let-7a-1:hsa-let-7a-3p", "hsa-let-7a-1:hsa-let-7a-5p", "hsa…
## $ FT42  <dbl> 444, 96079, 3, 96112, 450, 96864, 43, 7718, 3, 6968, 226, 137, 1…
## $ FT52  <dbl> 647, 164316, 2, 164372, 649, 165176, 78, 12420, 35, 14281, 476, …
## $ FT54  <dbl> 709, 293191, 22, 293225, 712, 294468, 294, 46712, 19, 39354, 828…
## $ FT60  <dbl> 813, 226991, 12, 227042, 823, 228053, 127, 31511, 74, 12717, 981…
## $ FT67  <dbl> 495, 186434, 8, 186507, 501, 187442, 66, 12363, 18, 23337, 292, …
## $ FT73  <dbl> 314, 55365, 3, 55392, 314, 55763, 25, 4229, 8, 6943, 161, 134, 1…
## $ FT85  <dbl> 276, 63794, 0, 63813, 278, 64435, 37, 7375, 7, 3376, 161, 99, 11…
## $ FT65  <dbl> 1493, 522367, 16, 522475, 1506, 524275, 214, 49933, 38, 20605, 1…
## $ FT82  <dbl> 1183, 396173, 3, 396201, 1196, 398205, 80, 25687, 33, 17536, 583…
## $ FT74  <dbl> 1483, 306320, 7, 306488, 1498, 308350, 160, 27529, 50, 13747, 98…
## $ PL166 <dbl> 37, 4892, 2, 4885, 37, 4904, 1, 573, 0, 763, 24, 49, 15, 741, 1,…
## $ PL165 <dbl> 377, 34289, 2, 34249, 380, 34388, 10, 4663, 13, 5342, 199, 272, …
## $ PL167 <dbl> 26, 4281, 0, 4278, 26, 4291, 3, 1333, 3, 972, 37, 65, 15, 1172, …
## $ FT36  <dbl> 82, 18877, 1, 18851, 82, 18994, 8, 2981, 3, 2978, 73, 169, 18, 3…
## $ FT41  <dbl> 161, 54184, 1, 54159, 163, 54412, 10, 5212, 0, 7149, 81, 226, 40…
## $ PM92  <dbl> 205, 53078, 5, 53113, 205, 53148, 49, 6631, 6, 3194, 167, 191, 6…
## $ PL149 <dbl> 373, 97396, 2, 97363, 374, 97813, 40, 6804, 3, 6575, 267, 364, 1…
## $ PL170 <dbl> 63, 8462, 1, 8458, 63, 8505, 1, 643, 5, 568, 25, 40, 27, 1884, 0…
## $ FT91  <dbl> 284, 36696, 0, 36658, 287, 36976, 21, 2567, 4, 2243, 190, 175, 9…
## $ PM84  <dbl> 801, 171395, 3, 171415, 803, 171989, 176, 39280, 18, 11563, 1017…
## $ PL148 <dbl> 717, 154473, 9, 154346, 723, 155008, 51, 8912, 18, 9547, 351, 50…
## $ PL171 <dbl> 43, 3864, 3, 3867, 44, 3871, 2, 524, 3, 262, 22, 24, 20, 1374, 1…
## $ PM114 <dbl> 1761, 382620, 5, 382570, 1765, 383478, 202, 39757, 48, 24359, 16…
## $ PM120 <dbl> 146, 31237, 0, 31180, 146, 31221, 20, 2763, 6, 1372, 70, 144, 55…
## $ PM167 <dbl> 737, 186118, 8, 185711, 745, 185863, 74, 13911, 15, 8704, 412, 9…
## $ PM181 <dbl> 120, 27350, 0, 27281, 122, 27375, 7, 2578, 2, 1433, 77, 158, 50,…
## $ PM135 <dbl> 845, 161419, 5, 161003, 851, 161174, 54, 14197, 22, 9396, 557, 6…
## $ PM153 <dbl> 613, 109686, 13, 109434, 614, 109331, 70, 13280, 21, 6478, 568, …
## $ PM112 <dbl> 1103, 214843, 4, 214606, 1104, 215051, 201, 29726, 26, 11124, 95…
## $ NTD4  <dbl> 944, 186632, 4, 186168, 959, 187011, 110, 12151, 29, 10555, 490,…
```

**CHECKPOINT:** Your filtered df, `eFil` should have 2344 rows and 31 columns  

#### **3.4 Remove duplicate mature miRNAs**     

Each mature miRNA (the ones having -3p and -5p at the end of their names) comes from precursor miRNAs - one precursor (e.g., miR-200) gives rise to two mature miRNAs (miR-200-3p and miR-200-5p). Now, there are a few precursors which give rise to the same two mature miRNAs (usually because it's only a few nucleotide difference somwehere in the hairpin loop of the precursor miRNA sequence).  

Let's first check which mature miRNAs are duplicated:  

**Separate the mirs column, count the number of mature miRNAs, and use filter to keep ones with a count of more than 1**


```r
eFil %>% 
    separate("mirs", c("precursor", "miRNA"), ":")  %>% 
    add_count (miRNA)  %>% 
    filter (n>1) 
```


So, we can see that there are 149 mature miRNAs that have been duplicated - We're going to keep the ones which have the highest overall expression  



```r
eFil_last <- eFil %>% 
  separate("mirs", c("precursor", "miRNA"), sep = ":") %>% 

  mutate(Sum = rowSums(select_if(., is.numeric))) %>% 
  # here I'm specifying to select only the numeric columns, because our first column are the miRNA names
  
  group_by(miRNA) %>% 
  # grouping by the mature miRNA names
  
  slice_max(Sum) %>% 
  # within the miRNA groups, only keeping the ones which have the highest value in Sum
  
  dplyr::select(-Sum, -precursor) %>%
  # removing the Sum and miRNA precursor column (as we have other dfs we can match to find out the precursor if needed)
  
  distinct(miRNA, .keep_all = TRUE)  %>% 
  glimpse()
```

```
## Rows: 2,137
## Columns: 31
## Groups: miRNA [2,137]
## $ miRNA <chr> "hsa-let-7a-2-3p", "hsa-let-7a-3p", "hsa-let-7a-5p", "hsa-let-7b…
## $ FT42  <dbl> 3, 450, 96864, 43, 7718, 3, 6968, 226, 137, 186, 11824, 3, 3, 32…
## $ FT52  <dbl> 2, 649, 165176, 78, 12420, 35, 14281, 476, 271, 316, 19141, 17, …
## $ FT54  <dbl> 22, 712, 294468, 294, 46712, 19, 39354, 828, 931, 194, 29262, 56…
## $ FT60  <dbl> 12, 823, 228053, 127, 31511, 74, 12717, 981, 605, 345, 25853, 13…
## $ FT67  <dbl> 8, 501, 187442, 66, 12363, 18, 23337, 292, 287, 266, 22143, 6, 8…
## $ FT73  <dbl> 3, 314, 55763, 25, 4229, 8, 6943, 161, 134, 115, 6131, 4, 7, 297…
## $ FT85  <dbl> 0, 278, 64435, 37, 7375, 7, 3376, 161, 99, 110, 14583, 0, 11, 44…
## $ FT65  <dbl> 16, 1506, 524275, 214, 49933, 38, 20605, 1805, 948, 464, 36042, …
## $ FT82  <dbl> 3, 1196, 398205, 80, 25687, 33, 17536, 583, 978, 339, 47495, 6, …
## $ FT74  <dbl> 7, 1498, 308350, 160, 27529, 50, 13747, 981, 665, 671, 36465, 16…
## $ PL166 <dbl> 2, 37, 4904, 1, 573, 0, 763, 24, 49, 15, 741, 1, 0, 447, 2, 219,…
## $ PL165 <dbl> 2, 380, 34388, 10, 4663, 13, 5342, 199, 272, 91, 7015, 4, 12, 51…
## $ PL167 <dbl> 0, 26, 4291, 3, 1333, 3, 972, 37, 65, 15, 1172, 0, 3, 1078, 11, …
## $ FT36  <dbl> 1, 82, 18994, 8, 2981, 3, 2978, 73, 169, 18, 3031, 2, 3, 3281, 1…
## $ FT41  <dbl> 1, 163, 54412, 10, 5212, 0, 7149, 81, 226, 40, 6757, 2, 3, 3399,…
## $ PM92  <dbl> 5, 205, 53148, 49, 6631, 6, 3194, 167, 191, 67, 6714, 4, 8, 6797…
## $ PL149 <dbl> 2, 374, 97813, 40, 6804, 3, 6575, 267, 364, 165, 14810, 2, 13, 6…
## $ PL170 <dbl> 1, 63, 8505, 1, 643, 5, 568, 25, 40, 27, 1884, 0, 2, 778, 14, 33…
## $ FT91  <dbl> 0, 287, 36976, 21, 2567, 4, 2243, 190, 175, 91, 7093, 4, 5, 3684…
## $ PM84  <dbl> 3, 803, 171989, 176, 39280, 18, 11563, 1017, 1004, 333, 34799, 1…
## $ PL148 <dbl> 9, 723, 155008, 51, 8912, 18, 9547, 351, 501, 305, 30924, 6, 17,…
## $ PL171 <dbl> 3, 44, 3871, 2, 524, 3, 262, 22, 24, 20, 1374, 1, 2, 698, 5, 244…
## $ PM114 <dbl> 5, 1765, 383478, 202, 39757, 48, 24359, 1672, 2063, 495, 47170, …
## $ PM120 <dbl> 0, 146, 31221, 20, 2763, 6, 1372, 70, 144, 55, 5139, 1, 13, 2536…
## $ PM167 <dbl> 8, 745, 185863, 74, 13911, 15, 8704, 412, 932, 259, 20695, 11, 1…
## $ PM181 <dbl> 0, 122, 27375, 7, 2578, 2, 1433, 77, 158, 50, 4232, 3, 2, 2357, …
## $ PM135 <dbl> 5, 851, 161174, 54, 14197, 22, 9396, 557, 688, 450, 22656, 5, 20…
## $ PM153 <dbl> 13, 614, 109331, 70, 13280, 21, 6478, 568, 802, 259, 13170, 3, 2…
## $ PM112 <dbl> 4, 1104, 215051, 201, 29726, 26, 11124, 950, 1281, 399, 27136, 1…
## $ NTD4  <dbl> 4, 959, 187011, 110, 12151, 29, 10555, 490, 821, 401, 24939, 14,…
```

```r
  # remove rows with duplicates within the miRNA column, but keep all the rest of the columns
  
  # column_to_rownames("miRNA")  # this decreased the number of the columns so I commented it out and I needed to keep them as columns for the next step too
```

**CHECKPOINT:** Your filtered df with mature miRNAs `eFil` should have 2137 rows and 31 columns 

**Make the density and boxplots for our filtered data, as we did for the raw data above**  


```r
eFil_plot <-  eFil_last %>% 
  pivot_longer(cols = -c("miRNA"),names_to = "Case_ID", values_to = "Reads") %>%  #I am not sure if I was supposed to do this step
  inner_join(., pDat, by = "Case_ID")

d2 <- eFil_plot %>% 
  ggplot(aes(x = log10(Reads + 1), colour = Case_ID)) +
  geom_density() +
  theme_minimal() +
  theme(legend.position = "none")

b2 <- eFil_plot %>% 
  ggplot(aes(x = Case_ID, y = log10(Reads + 1), fill = Case_ID)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none")  +
  labs(title = "Gene Expression", subtitle = "Raw Counts") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggarrange(d1, b1, d2, b2, nrow = 2, ncol = 2)
```

![](2_QC_files/figure-html/fil-plots-1.png)<!-- -->

***

### **4.0 Normalization**  

Now that we've filtered our data, it's not in the 'shape' it was before - you can also see for the plots above that the maximum density is also low (which makes sense, and is good confirmation since we removed miRNAs with 0 exp in all samples).  

To take into account that the relative expression (exp between one miRNA to another) between the miRNAs within our eFil dataset is different from that of the relative expression between miRNas in our original, raw `eDat` df, we use a techniques known as normalization  

Normalization takes into account the full expression dataset we have, and adjusts the expression of each miRNA, realtive to the whole dataset. It outputs a numerical value for each miRNA which is either subtracted to and or added to the original expression  

We're going to use the [Relative Log Expression](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0206312#:~:text=Relative%20Log%20Expression%20(RLE)%3A,geometric%20mean%20across%20all%20samples.) method of normalization


```r
eFil <- eFil %>% 
column_to_rownames(var="mirs")

filteres_mir_names <- as.data.frame(row.names(eFil))
samples <- pDat$Case_ID

# computing normalized counts (between sample normalization of miRNA expression)
norm_counts <- DGEList(counts = eFil, samples = samples, genes = filteres_mir_names)
norm_counts <- calcNormFactors(norm_counts, method = "RLE")

# performing within sample normalization of miRNA expression (CPM = Counts per Million = RPM = Reads per Million)
eNorm <- cpm(norm_counts)

# converting to our usual log10 expression
eNorm <- log10(eNorm + 1) 

# converting the expression counts into whole numbers
eNorm <- round(eNorm, digits = 0)

eNorm <- as.data.frame(eNorm)
```

**Make the density and boxplots for our normalized data, as we did for the raw data above**  


```r
eNorm_plot <- eNorm %>% 
  rownames_to_column(var="miRNA") %>% 
  pivot_longer(cols = -c("miRNA"),names_to = "Case_ID", values_to = "Reads")
  

d3 <- eNorm_plot %>% 
  ggplot(aes(x=Reads, colour = Case_ID)) +
  geom_density() +
  my_theme +
  theme(legend.position = "none")

b3 <- eNorm_plot %>% 
  ggplot(aes(x = Case_ID, y = Reads, fill = Case_ID)) +
  geom_boxplot() +
  my_theme +
  theme(legend.position = "none")  +
  labs(title = "Normalized Gene Expression")

ggarrange(d1, b1, d2, b2, d3, b3, nrow = 3, ncol = 2)
```

![](2_QC_files/figure-html/norm-plots-1.png)<!-- -->

As we can see now, the expression within our normalized df `eNorm` is much well distributed than our raw expression  

***

### **5.0 Sample-sample correlation**  

As an additional measure, we'll check how well our samples relate to each other - i.e. are samples within the same trimester showing higher degrees of association (which is what we would expect)


```r
raw_cor <- cor(eDat)

norm_cor <- cor(eNorm)

# for the following function, we require the column names of our expression dataframe to be the rownames of our pDat dataframe. Otherwise the function doesn't work
pDat <- pDat %>% 
  column_to_rownames("Case_ID")


raw_cor %>% 
  pheatmap::pheatmap(clustering_distance_cols = "euclidean", clustering_method = "complete", 
                     #the method by which we're caluclating the `distance` i.e. the similarity between the samples. The `complete` method of clutering deals better with outlier samples than outer methods
                     
                     cluster_rows = TRUE, 
                     #clustering by rows
                     
                     show_colnames = TRUE, show_rownames = TRUE,
                     #show sample names for both rows and columns
                     
                     annotation_col = pDat[c("Trimester")], annotation_row = pDat[c("Trimester")],
                     #select the variables by which we want to observe our clustering
                     
                     #name of the plot
                     main = "Raw Sample Correaltion") 
```

![](2_QC_files/figure-html/cor-maps-1.png)<!-- -->

```r
norm_cor %>% 
  pheatmap::pheatmap(clustering_distance_cols = "euclidean", clustering_method = "complete", 
                     cluster_rows = TRUE, show_colnames = TRUE, show_rownames = TRUE,
                     annotation_col = pDat[c("Trimester")], annotation_row = pDat[c("Trimester")],
                     main = "Raw Sample Correaltion")
```

![](2_QC_files/figure-html/cor-maps-2.png)<!-- -->

This shows us that by normalizing, samples within the same trimester cluster together and have high correlation between each other  

***  

### **5.0 Dimensionality Reduction**  

Dimensionality reduction techniques allow us to observe 'big data' in a meaningful way, where all of the expression data per sample and within sample is converted into units showing the 'distance' i.e., similarity between them. Here's we'll be using Principal Component Analysis, which reduces and converts our data points into distinct Principal Components (which have no units of measure) and shows us how much variability exists between our data, and which variables are contributing towards that variation  


```r
# first convert eNorm to log10+1
eLog <- log10(eNorm + 1)

# the function we use performs it's operation row-wise, and since we want to check the variability between genes and not the samples, we'll transpose our exp df to convert the samples to be rows, and the columns to be our genes. Don't open this df as it might cause R to crash (as there will be more columns than rows in this df now)
t_eLog <- as.data.frame(t(eLog))

# the function to compute PCA. We set scale to FALSE, because we already log10 scaled our data, and center to TRUE, so that the values are converted to Z-scores
pca <- prcomp(t_eLog, scale = FALSE, center = TRUE) #n has to be less than number of samples

# each score represents the contribution of each gene within each sample
scores <- pca$x

summary(pca)
```

```
## Importance of components:
##                           PC1    PC2     PC3     PC4     PC5     PC6     PC7
## Standard deviation     0.9687 0.9479 0.67465 0.55795 0.54592 0.49771 0.47483
## Proportion of Variance 0.1514 0.1450 0.07344 0.05023 0.04809 0.03997 0.03638
## Cumulative Proportion  0.1514 0.2964 0.36982 0.42006 0.46815 0.50812 0.54450
##                            PC8     PC9    PC10    PC11    PC12    PC13    PC14
## Standard deviation     0.46861 0.46570 0.43848 0.41899 0.41672 0.40033 0.38815
## Proportion of Variance 0.03543 0.03499 0.03102 0.02833 0.02802 0.02586 0.02431
## Cumulative Proportion  0.57993 0.61493 0.64595 0.67428 0.70230 0.72816 0.75247
##                          PC15    PC16    PC17    PC18    PC19   PC20   PC21
## Standard deviation     0.3784 0.36854 0.36203 0.34441 0.33721 0.3274 0.3217
## Proportion of Variance 0.0231 0.02192 0.02115 0.01914 0.01835 0.0173 0.0167
## Cumulative Proportion  0.7756 0.79749 0.81864 0.83778 0.85613 0.8734 0.8901
##                           PC22    PC23    PC24   PC25    PC26    PC27    PC28
## Standard deviation     0.31473 0.30800 0.30280 0.3008 0.29150 0.28003 0.27200
## Proportion of Variance 0.01598 0.01531 0.01479 0.0146 0.01371 0.01265 0.01194
## Cumulative Proportion  0.90611 0.92142 0.93621 0.9508 0.96452 0.97717 0.98911
##                           PC29      PC30
## Standard deviation     0.25977 8.159e-16
## Proportion of Variance 0.01089 0.000e+00
## Cumulative Proportion  1.00000 1.000e+00
```

Here, the summary output shows us the how much variability each PC is contributing in the Proportion of Variance row, and the total sum variability each of each preceding PC in the Cumulative Proportion row. The Cumulative Proportion will always add up to a 100% at the last PC.  

> We're now going to make the 3rd type of file you'll usually find along with `eDat` and `pDat` - the metadata, or `mDat` file.  

This file contains additional information about our samples, that isn't clinical/phenotypical in nature  


```r
scores <- scores %>% 
  as.data.frame() %>% 
  rownames_to_column("Case_ID")

pDat <- pDat %>% 
  rownames_to_column("Case_ID")

mDat <- pDat %>% 
  left_join(scores, by = "Case_ID")



# getting the data from the summary function we saw above  
summary <- data.frame(PC = 1:30, 
                      var_explained = (pca$sdev)^2 / sum((pca$sdev)^2), 
                      cumulative = cumsum(pca$sdev^2 / sum(pca$sdev^2))
                      )

summary <- summary %>% 
  mutate(cumulative_perc = cumulative*100)

# we're now going to convert our PC column into a different data structure - a factor. Factors allow you to order variables within a column according to your liking; each variable is now called a level

# the sort argument is telling to order the column by ascending numerical order
summary <- summary %>% 
  mutate(PC = sort(as.factor(PC)))

summary$PC
```

```
##  [1] 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
## [26] 26 27 28 29 30
## 30 Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 ... 30
```

**Convert `Trimester` in pDat to a factor**  
*Hint: Use the `fct_relevel` function*


```r
pDat <- pDat %>% 
  mutate (Trimester = fct_relevel(Trimester))
```

**Make plots with PCs on the x-axis and** 

- var_explained on the y-axis  
- cumulative_perc on the y-axis


```r
summary %>% 
  ggplot(aes(x = PC, y = var_explained)) +
  geom_boxplot() +
  my_theme +
  theme(legend.position = "none")  +
  labs(title = "Normalized Gene Expression")
```

![](2_QC_files/figure-html/pc-plots-1.png)<!-- -->

We see from the summary and the plots, that PC1 is contributing about 15% of the variance, PC2 13%, and so on, to add up to 100% at PC30 (30 PCs because we have 30 samples)  

We'll now check which of our variables actually contribute towards the PCs and drive variance within each PC.  

But first, we'll select the variables of interest: Trimester, Sex, Processing_time, RQS, flowCell   

Why flowCell? My previous analyses have shown that [flowcell](https://www.youtube.com/watch?v=pfZp5Vgsbw0)(which is a kind of chip on which sequencing is carried out. Each flowcell has up to 8 lanes, and a sample is run per lane) contributes towards the variation we see, and we get different differentially-expressed genes depending upon if we include it or not. To quickly check of it indeed does, I'll use as an example  


```r
eDat_plot %>% 
  separate("miRNA", c("precursor", "miRNA"), ":") %>%   #not sure why at this point the column wasn't seperated into two so had to add 
  filter(miRNA == "hsa-miR-21-5p") %>% 
  ggplot(aes(x = miRNA, y = log2(Reads), fill = Trimester)) +
  geom_boxplot() +
  facet_grid(~flowCell) +
  theme_minimal() 
```

![](2_QC_files/figure-html/flowcell-1.png)<!-- -->

Here, we can see that on two of the flowcells, only trimester 2 samples were run, whereas on the last flow cell, we have a mix of all trimesters. As you can imagine, not having the samples equally distributed amongst all 3 flow cells means that cells 1 and 2 are biased towards expression for trimester 2 samples.  


```r
# Selecting the variables of interest
pDat_cor <- pDat %>% 
  dplyr::select(Trimester, Sex, Processing_time, RQS, flowCell)

scores <- scores %>% 
  column_to_rownames("Case_ID")

# calculating the pvalue significance of the above variables by PC
pc_pval <- plomics::lmmatrix(dep = scores, ind = pDat_cor, metric = "Pvalue") #correaltion matrix

# transposing
pc_pval <- t(pc_pval) #only for visual purposes - Not being used for any calculations

# highlighting the variables significantly contribuing towards each PC
pc_pval %>% 
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
   <th style="text-align:left;"> Processing_time </th>
   <th style="text-align:left;"> RQS </th>
   <th style="text-align:left;"> flowCell </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.579</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.001</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.12</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.707</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.39</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.891</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.852</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.376</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.001</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.276</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.002</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.936</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.223</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.661</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.222</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.019</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.65</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.33</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.934</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.527</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.596</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.806</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.454</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.789</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.519</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.77</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.834</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.833</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.411</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.295</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.789</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.886</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.131</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.487</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.127</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.634</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.832</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.887</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.355</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.382</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.142</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.985</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.768</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.957</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.712</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.846</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC11 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.979</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.384</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.956</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.826</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">1</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC12 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.982</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.725</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.792</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.285</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.178</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC13 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.759</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.311</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.47</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.828</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.644</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC14 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.957</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.04</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.3</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.282</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.813</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC15 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.813</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.485</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.431</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.426</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.774</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC16 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.958</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.466</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.106</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.551</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.794</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC17 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.91</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.329</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.664</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.868</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.636</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC18 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.955</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.418</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.754</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.802</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.853</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC19 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.919</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.709</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.564</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.734</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.508</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC20 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.984</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.884</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.925</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.891</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.318</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC21 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.989</span> </td>
   <td style="text-align:left;"> <span style="     color: green !important;">0.043</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.497</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.71</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.93</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC22 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.921</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.835</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.961</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.749</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.995</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC23 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.961</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.092</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.716</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.865</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.995</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC24 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.965</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.636</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.958</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.96</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.979</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC25 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.946</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.215</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.937</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.427</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.901</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC26 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.978</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.54</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.642</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.908</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.055</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC27 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.98</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.066</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.812</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.926</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.963</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC28 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.957</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.647</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.99</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.49</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.987</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC29 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.941</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.511</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.827</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.996</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.523</span> </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC30 </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.454</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.059</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.775</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.471</span> </td>
   <td style="text-align:left;"> <span style="     color:  !important;">0.1</span> </td>
  </tr>
</tbody>
</table>


We're now going to represent the data in the above table in a plot:  



```r
pc_plot <- pc_pval %>% 
  as.data.frame %>% 
  mutate(Principal_Component = rownames(pc_pval),
         PC = as.factor(1:30)) %>% 
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

pc_plot %>% 
  ggplot(aes(x = PC, y = Variable , fill = pval_cat)) +
  geom_tile(col = "lightgray") +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + #expand function fits the plot to its assigned dimesions 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = pc_colpal)  + 
  coord_fixed() + # very important otherwise height of each tile is elongated!
  labs(title = "Contribution of Variables to Observed Expression\n", x = "", y = "", fill = "P value")
```

![](2_QC_files/figure-html/pval-plot-1.png)<!-- -->

We'll do the same to calculate the R-saured value for each of the variables within each PC.  

**Can you explain what an R-squared value is?**  
It is a value that showed the degree of correlation.it ranges between 0 and 1. closer to 1, stronger is the correlation. 


```r
pc_rsq <- lmmatrix(dep = scores, ind = pDat_cor, metric = "Rsquared") #correaltion matrix

pc_rsq <- t(pc_rsq)

rsq_plot <- pc_rsq %>% 
  as.data.frame %>% 
  mutate(Principal_Component = rownames(pc_rsq),
         PC = as.factor(1:30)) %>% 
  pivot_longer(cols = -c(Principal_Component, PC), names_to = "Variable", values_to = "rsq")

rsq_plot %>% 
  ggplot(aes(x = PC, y = Variable , fill = rsq)) +
  geom_tile(col = "lightgray") +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + #expand function fits the plot to its assigned dimesions 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "#509e36") +
  coord_fixed() + # very important otherwise height of each tile is elongated!
  labs(x = "", y = "", fill = "Rsq")
```

![](2_QC_files/figure-html/rsq-plot-1.png)<!-- -->

Viewing the plots together:


```r
p1 <- pc_plot %>% 
  ggplot(aes(x = PC, y = Variable , fill = pval_cat)) +
  geom_tile(col = "lightgray") +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + #expand function fits the plot to its assigned dimesions 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = pc_colpal)  + 
  coord_fixed() + # very important otherwise height of each tile is elongated!
  labs(title = "Contribution of Variables to Observed Expression\n", x = "", y = "", fill = "P value")

r1 <- rsq_plot %>% 
  ggplot(aes(x = PC, y = Variable , fill = rsq)) +
  geom_tile(col = "lightgray") +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + #expand function fits the plot to its assigned dimesions 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "#509e36") +
  coord_fixed() + # very important otherwise height of each tile is elongated!
  labs(x = "", y = "", fill = "Rsq")

egg::ggarrange(p1, r1, nrow = 2, ncol = 1, heights = c(3,3)) 
```

![](2_QC_files/figure-html/p-rsq-1.png)<!-- -->

Now, we can confidently say that of the 5 variables we selected, barring Sex, all the other 4 are contributing towards the variation in expression that we see, and so we'll be including them further on in our analyses  

Now that we've figured which variables to use, let's actually check how the samples separate by them  

**Make the following scatter plots with mDat**  

- PC1 vs PC2 for:
  - Trimester
  - RQS (bucket according to your preference)
  - hours
  - flowcell

- PC2 vs PC3 for:
  - Trimester
  - RQS
  
Colour by variables mentioned above, and set the same x and y axes +/- limits  


```r
pc1_pc2_1 <- mDat %>% 
  ggplot(aes(x = PC1, y = PC2 , color = Trimester)) +
  geom_point()+
  scale_x_continuous(limits = c(-2, 2)) + 
  scale_y_continuous(limits = c(-2, 2))
  

pc1_pc2_2 <- mDat %>% 
  ggplot(aes(x = PC1, y = PC2 , color = RQS)) +
  geom_point()+
  scale_x_continuous(limits = c(-2, 2)) + 
  scale_y_continuous(limits = c(-2, 2))
  
pc1_pc2_3 <- mDat %>% 
  ggplot(aes(x = PC1, y = PC2 , color = Processing_time)) +
  geom_point()+
  scale_x_continuous(limits = c(-2, 2)) + 
  scale_y_continuous(limits = c(-2, 2))


pc1_pc2_4<- mDat %>% 
  ggplot(aes(x = PC1, y = PC2 , color = flowCell)) +
  geom_point()+
  scale_x_continuous(limits = c(-2, 2)) + 
  scale_y_continuous(limits = c(-2, 2))

ggarrange(pc1_pc2_1, pc1_pc2_2,pc1_pc2_3, pc1_pc2_4, nrow = 2, ncol = 2)
```

![](2_QC_files/figure-html/PC1-2-1.png)<!-- -->



```r
pc2_pc3_1 <- mDat %>% 
  ggplot(aes(x = PC2, y = PC3 , color = Trimester)) +
  geom_point()+
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2))

pc2_pc3_2 <- mDat %>% 
  ggplot(aes(x = PC2, y = PC3 , color = RQS)) +
  geom_point()+
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2))

ggarrange(pc2_pc3_1, pc2_pc3_2, nrow = 2, ncol = 1)
```

![](2_QC_files/figure-html/PC2-3-1.png)<!-- -->

Based on the plots, which variable/s seem to be contributing the most variation towards our data?  

***


