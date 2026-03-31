
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

Simulation of proteomic data matrices were used to test the behavior of each 
item in the assessment of metric and graphical information, and also to
eval the global score. 

Function structure:

* A first block of code corresponds to the generation of plain proteomic data. Two
groups of study are considered, and different parameters can be changed, including 
number of proteins; internal correlation between samples and proteins; proportion 
of proteins with different quantities between groups; residual error; and
quantification range.

* With the aim of increasing variability in data to mimic the different behavior 
that must be detected across the score items, two extra blocks of code allow the 
addition of a shift and a dependence between mean and intensity in the data. 

* Finally, the last block of code generates missing data following a pattern of 
low quantity (not random missing data). 

Using the predefined parameter values, the function simulates data with perfect 
characteristics for all items. By changing parameters, error
is included into the data, with the corresponding behaviour change. 

For each item, nine parameter values were selected, which was called "specific
combination". The first value was used to 
generate a perfect dataset, the ninth to generate the worst dataset, and the 
remaining values to generate intermediate datasets with progressively worse quality.
The item was then used to rate and sort the datasets, and Kendall correlation
was computed between the observed and expected order. Nine specific combinations
were stablish to recreate and assess different situations regarding datasets 
behaviour. 

For each specific combination, 9 dataset were simulated and compared using the item, 
obtaining one correlation. Then, these simulations were repeated for a different 
combination of seeds (10), number of proteins (4) and sample size per group (4). 
Therefore, a total number of 1440 combinations were obtained:

$$\text{9 specific combinations} \times \text{10 repetitions} \times \text{4 diff nº of proteins}  \times \text{4 different sample size} = 1440$$

For each specific combination, 9 datasets were simulated: 12960 simulations. 

For general combination, one Kendall correlation is computed = 1440 Kendall correlation values. 














