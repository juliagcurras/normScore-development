
# normScore

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
The output is a correction factor with a value between 1 and 0. 

* For each data matrix normalized with a different method, items 1 to 6 are obtained. 
The input is the normalized data matrix, but sometimes group infomation is also required. 
The output is a different value for each item. 

* Once an item is estimated, the different normalizations are ranked based on the 
item values. The highest values corresponds to the poorest performance. 

* Finally, the sum of ranking positions for each normalization is computed. A correction
factor for taking into account the normalization requirement is applied to the 
logarithm method (non-normalised data). Then, the normalization with the lowest 
score is selected as the best one. 

In the following sections, the estimation of each item is explaining. 


## Item 0 - correction factor

To assess the magnitude of the systematic bias, total intensity (sum of protein 
intensities) is estimated for each sample. Then, the maximum values are
compared with the minimal ones. An strong difference between the previous values 
indicates the presence of important variability in the data, that has to be removed 
using a normalization method. 


In this line, the estimation of the correction factor is based on the magnitude 
of the systematic bias present in the data. Steps:

1. Using the sum of total intensity in each sample, the ratio between the lowest 
value and each of the rest intensities is calculated, getting a group of ratios.

2. Then, the median value of the previous ratios is estimated.

3. Step 1 and 2 are repeated for: A) the 3 lowest values of total intensity when 
the number of samples < 16; B) the 30% of minimal values when the number of samples
is equal to or higher than 16. 

4. The mean of the ratios generated in the previous steps (one ratio for each 
minimum value of total intensity considered in step 3) is calculated. 

5. The correction factor is obtained as 1 - mean from step 4. 


The score for the non-normalized data (only log transformed) is corrected using 
this factor, so:

* If the final ratio between minimal and other values is close to 1, the correction
factor will be close to 0, so after the correction, the score for non-normalized
data will decrease, climbing on the ranking positions. A ratio close to 1 means that
there is no big differences between total intensities, which reflecs a 
lack of systematic biased. 

* If the final ratio is close to 0, the correction factor will be close to 1, so
after the correction the score for non-normalized data will remain practically 
the same. A final ratio similar to 0 reflects an important systematic bias in the 
data, so applying a normalization method to remove that variations is necessary. 

A final ratio is 0.5 or lower means that minimum intentisities are half or 
more than half of the maximum intensities. In this situation, the corrected factor 
become 1 automatically, because a normalization is needed and there is no reason
to decrease the score of only log-transformed data. 



## Item 1

The coefficiente of variance (CV) is estimated for each sample across their protein 
intensity values, and then the mean of the CV is obtained for each group of samples.
Then, the mean between groups is estimated, being this value the one used to rank
the different normalizations in item 1. 

## Item 2

For each pair of samples from the same group, Spearman correlation is estimated 
between protein quantities. Then, the median for the different correlation is
estimated and is corrected usig a third part of the IQR (interquartilic range) 
value. To meet the requirement: higher values, lower performance, this final value
is substracted from 1, and this is the score of item 2. 

## Item 3

For item 3, logFC and sum of average intensity in each of the defined groups are
estimated and then used to adjust a regression line. Ideally, this line should 
have intecerpt = 0, coefficiente of regression = 0. The difference between the 
expected line and the observed line is estimated, obtaining and area. This area 
needs a standarization, that is done using the min-max values (previously used 
to estimated the difference between regression lines). The standarized area is 
the value used as item 3 to rank the normalizations. 

## Item 4

Similar to item 3, in these case there should be an independence between standard
deviation and samples sorted by the mean of intensities of each sample. Therefore, 
a regression line is estimated using the previous data, and then is compared with
the expected line (intecerpt = 0, coefficient of regression = 0), The area is also
standarized using the min-max values, and the final standarized area is the value
for item 4. 

## Item 5

If we estimate the ratio of each protein quantity and the mean of quantity in each 
sample) and then we transform these values with the logarithm, we expected a median
for each sample of 0. Why? Because the correspondent value without transformation 
is 1 which means that ratio protein quantity = average quantity across samples. 

Thus, we expected a median value in each of the samples for this log-transformed 
ratio of 0. So we estimate the mean squared error (MSE) between observed and expected
medians, and the normalization we the lowest MSE will be the best one. 

## Item 6

The distribution of the protein quantities across samples should be aligned in 
the perfect normalization. Thus, the following mse are computed:

* Median of protein quantity in each samples (observed) vs global median (expected)

* Q1 of protein quantity in each samples (observed) vs global Q1 (expected)

* Q3 of protein quantity in each samples (observed) vs Q3 median (expected)
























