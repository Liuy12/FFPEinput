# FFPEinput

The tutorial aims to provide RNA input recommendations (RNA concentration and library concentration) for FFPE samples in order to generate adequate data for RNA-seq. Please refer to our [publication](https://assets.researchsquare.com/files/rs-1428294/v1/5f199ab7-8b30-4cf0-a3eb-296477271c8b.pdf?c=1647869411) for more details

## Background

Formalin-fixed, paraffin-embedded (FFPE) tissues have many advantages for identification of risk biomarkers, including wide availability and potential for extended follow-up endpoints. However, RNA derived from archival FFPE samples has limited quality. There are limited studies that provide insight for selection of FFPE samples of adequate quality. Here we identified parameters that determine which FFPE samples have the potential for successful RNA extraction, library preparation, and generation of usable RNAseq data.

## Installation

All analysis is performed in R. If you don't have it, you can download R from [here](https://www.r-project.org/). Then you will need to install required R package.

```r
install.packages("caret")
```

## Data

You can download the .txt file within the data folder under this repo. Then load it within a R session.

```r
sampleinfo <- read.delim('path/to/file',header=T,stringAsFactors=F,check.names=F)
str(sampleinfo)
```

The data contains pre-sequencing lab metrics and post-sequencing bioinformatics metrics for 101 samples: 

- RNA ng/ul: input RNA concentration 
- Library ng/ul: pre-capture library concentration
- Gene mapped reads: total number of reads mapped to genic regions
- TPM > 4: total number of genes with Counts Per Million (TPM) bigger than 4. You can see the definition of TPM [here](https://www.reneshbedre.com/blog/expression_units.html)
- Median cor: spearman median correlation with the rest of samples within the cohort

## Define QC pass or fail based on bioinformatics metrics

QC status is determined based on the three bioinformatics metrics (Gene mapped reads, CPM > 2 & Median cor). 

```r
sampleinfo$FILTER <- ifelse(sampleinfo$`Gene mapped reads` > 25*10^6 & sampleinfo$`TPM > 4` > 11,400 & sampleinfo$`Median cor` > 0.75,'PASS','FAIL')
```

User are encouraged to choose cutoffs based on their preference of stringency.

## Bulid decision tree model using pre-sequencing lab metrics

```r
sed.seed(1234)
sample_split <- sample.split(Y = sampleinfo$FILTER, SplitRatio = 0.7)
train_set <- subset(x = sampleinfo, sample_split == TRUE)
test_set <- subset(x = sampleinfo, sample_split == FALSE)
modelfit <- train(FILTER ~ `RNA ng/ul` + `Library ng/ul`, method='rpart',data=train_set,
                   trControl = trainControl(method = "repeatedcv",repeats = 3),tuneLength=10)
plot(modelfit)
rpart.plot::rpart.plot(modelfit$finalModel, extra = 106)
```
