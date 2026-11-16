
# normScore

## Description

Score development for assessing normalization in proteomic/transcriptomic data.

The score consist of 6 different items that represents: 

  2) Item 1-2: a summary of pooled metrics for evaluating variability between samples
  
  1) Item 3-6: a numeric assessment of graphical representations; 
  
An individual score is estimated for each of the applied normalization methods. 
Then, a correction factor is applied to reflect the magnitude of systematic bias, 
that is, the necessity of normalization. 

The normalization with the lowest scores will be the one with the best performance. 

Main features:

* The input for item 0 estimation is the raw data matrix, without normalization or transformation. 
The output is a correction factor with a value from 0 to infinite. WARNING! Data must be previously
filtered by missing data. Missing data seriously affect this parameter. 

* For each data matrix normalized with a different method, items 1 to 6 are obtained. 
The input is the normalized data matrix, but sometimes group information is also required. 
The output is a different value for each item. 

* Once estimated for all normalization, the six items are scaled to 0-1 using the
min - max normalization, so in each column values will range from 0 to 1. 

* The sum of scaled values for each normalization is computed. The highest 
values corresponds to the poorest performance. A correction factor for taking into
account the normalization requirement is applied to the logarithm method 
(non-normalised data). 

* Finally, a 95% cofident interval is estimated for each score using a bootstrap
resample strategy (number of resampling = 1000) and a forest plot is created using
the scores and their interval. 


The simplest normalization with the lowest score is selected as the best one. 

In the following sections, the estimation of each item is explaining. 


### Item 0 - correction factor

To assess the magnitude of the systematic bias, total intensity (sum of protein 
intensities) is estimated for each sample. Then, the coefficiente of variance is 
estimated for this data. This coefficiente is used as the correction factor.

### Item 1

The coefficiente of variance (CV) is estimated for each sample across their protein 
intensity values, and then the mean of the CV is obtained for each group of samples.
Then, the mean between groups is estimated, being this value the one used to rank
the different normalizations in item 1. 

### Item 2

For each pair of samples from the same group, Spearman correlation is estimated 
between protein quantities. Then, the median for the different correlation is
estimated and is corrected usig a third part of the IQR (interquartilic range) 
value. To meet the requirement: higher values, lower performance, this final value
is substracted from 1, and this is the score of item 2. 

### Item 3

For item 3, logFC and sum of average intensity in each of the defined groups are
estimated and then used to adjust a regression line. Ideally, this line should 
have intecerpt = 0, coefficiente of regression = 0. The difference between the 
expected line and the observed line is estimated, obtaining an area. This area 
needs a standarization, that is done using the min-max values (previously used 
to estimated the difference between regression lines). 

For taking into account 
the expected shape, values are divided into 10 groups with the same sample size 
by average intensity. Then, Spearman correlation between IQR for logFC
in each group and the expected order by magnitude (decreasing) are computed. 
Correlation values, standarized to 0-1, are then used to correct the areas, 
estimated previously.

This corrected, standarized area is the value used as item 3 to rank the normalizations. 

### Item 4

Similar to item 3, in these case there should be an independence between standard
deviation and samples sorted by the intensities mean of each sample. Therefore, 
a regression line is estimated using the previous data, and then is compared with
the expected line (intecerpt = 0, coefficient of regression = 0), The area is also
standarized using the min-max values, and the final standarized area is the value
for item 4. 

### Item 5

If we estimate the ratio of each protein quantity and the mean of quantity in each 
sampl, and then we transform these values with the logarithm, we expected a median
for each sample of 0. Why? Because the correspondent value without transformation 
is 1 which means that ratio protein quantity = average quantity across samples. 

Thus, we expected a median value in each of the samples for this log-transformed 
ratio of 0. So we estimate the mean absolute percentage error (MAPE) between observed 
and expected medians of non log-transformed data, and the normalization with the 
lowest MAPE will be the best one. 

### Item 6

The distribution of the protein quantities across samples should be aligned in 
the perfect normalization. Thus, the following MAPE are computed:

* Median of protein quantity in each samples (observed) vs global median (expected)

* Q1 of protein quantity in each samples (observed) vs global Q1 (expected)

* Q3 of protein quantity in each samples (observed) vs Q3 median (expected)




## Testing

### Simulations

Simulated proteomic datasets were used to evaluate the behavior of each scoring 
item, as well as the overall performance of the `normScore` framework.

The simulation function is structured in three main components:

* **Baseline data generation**
  A base proteomic dataset is generated assuming two biological groups. 
  Several parameters can be controlled, including:

  * number of proteins
  * sample size per group
  * correlation structure between samples and proteins
  * proportion of differentially abundant proteins
  * residual noise
  * quantification range

* **Systematic perturbations**
  Additional components are introduced to mimic specific types of distortions 
  targeted by the score:

  * global intensity shifts
  * mean–variance dependence

* **Missing data generation**
  Missing values are introduced following a low-abundance (MNAR-like) 
  mechanism rather than random missingness.

Using the default parameters, the function generates datasets with ideal 
characteristics for all score items. By progressively modifying parameters, 
controlled levels of error are introduced, resulting in datasets of decreasing 
quality.


#### Item-wise evaluation

For each scoring item, a set of nine parameter values was defined 
(referred to as a *specific combination*), representing a gradient from 
optimal to poor data quality:

* the first value generates a near-perfect dataset
* the ninth value generates a highly distorted dataset
* intermediate values produce progressively degraded datasets

Each item was then used to rank the nine simulated datasets, and the agreement 
between the expected and observed ranking was assessed using Kendall’s correlation.


#### Simulation design

For each specific combination:

* 9 datasets were generated
* 1 Kendall correlation value was computed

This process was repeated across multiple simulation conditions:

* 20 random seeds
* 4 different numbers of proteins
* 4 sample sizes per group
* 6 levels of error magnitude

This results in:

$$
6 \times 20 \times 4 \times 4 = 1920 \text{ specific combinations}
$$

Since each combination contains 9 datasets:

* **Total datasets generated:** 17280
* **Total Kendall correlations computed:** 1920


#### Error types evaluated

Different types of perturbations were used to assess the sensitivity of each item:

* Items 0, 1, 5, 6: global intensity shift
* Item 2: correlation structure
* Item 3: shift + residual variance
* Item 4: mean–SD dependence
* normScore: combination of all error types

Considering all error scenarios:

* **Total datasets simulated:** 86,400
* **Total specific combinations:** 9,600


> These simulations allowed us to verify that each item is sensitive to the 
specific type of distortion it is designed to capture, and that the combined 
score behaves consistently across heterogeneous scenarios.

### One-hundred datasets

To evaluate the performance of the `normScore` framework beyond simulations, a
curated gold-standard dataset was assembled. This dataset consisted of 100 
proteomics studies from PRIDE, for which eight normalization methods were applied.

Each dataset was manually annotated twice, by two independent reviewers, 
by ranking normalization methods according 
to each individual criterion, and by identifying the best-performing methods 
overall. In some cases, clearly suboptimal methods were also flagged for 
exclusion. Then, normScore was also applied to the collection to get the 
item and global level rankings. 

Results were compared between normScore and manual reviewers using Spearman's 
correlation and pairwise correlation for item-levels rankings, as well as the 
Jaccard index for global rankings. 


### Benchmarking datasets

Finally, we used three publicly available proteomics benchmark datasets with 
known abundance ratios to assess how closely normScore rankings agreed with 
rankings based on independent measures of normalization performance.

- **PXD028735** contains DDA and DIA data from mixtures of human, 
*Escherichia coli*, and yeast proteins. The expected B/A fold changes (FCs) 
were 1, 4, and 0.5, respectively. We reanalysed the raw DDA data with MaxQuant 
v2.6.4 and the raw DIA data with DIA-NN v2.7.
- **PXD001819** contains yeast samples spiked with human UPS1 proteins and 
quantified using MaxQuant. We retained comparisons with expected UPS1 FCs of 
2 and 10. Yeast protein abundances were expected to remain unchanged.
- **PXD039399** contains four groups of yeast samples with protein inputs 
of 200, 100, 50, and 25 ng, quantified using DIA-NN. The expected FCs for 
the 200-ng group relative to the 100-, 50-, and 25-ng groups were 2, 4, and 8, 
respectively.

We preprocessed each dataset separately and applied the normalization methods 
under evaluation. To assess performance against the known experimental ratios, 
we calculated the mean absolute error (MAE) between estimated and expected FCs. 
Where applicable, we also calculated the area under the receiver operating 
characteristic curve (AUC) to assess how well the analysis distinguished proteins 
expected to be differentially abundant from those expected to remain unchanged.
We then used Kendall’s rank correlation to compare the normalization method 
rankings derived from each benchmark metric with those obtained using normScore.

For **PXD039399**, we computed normScore rankings both across all samples and 
separately within each group. For each pairwise comparison, we averaged the
rankings from the two groups involved.


---