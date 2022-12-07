FAQ
========

**1. Can I run a GREML analysis in a small sample?**

It is not recommended to run a GCTA-GREML analysis in a small sample. When the sample size is small, the sampling variance (standard error squared) of the estimate is large (see GCTA-GREML power calculator), so the estimate of SNP-heritability (h2-SNP) will fluctuate a lot and could even hit the boundary (0 or 1). Therefore, when the sample size is small, it is not surprising to observe an estimate of SNP-heritability being 0 or 1 (with a large standard error).

If the estimate hits the boundary (0 or 1), the phenotypic variance-covariance matrix (V) will often become invertible and you will see error message

.. code:: sh

    "Error: the variance-covaraince matrix V is not positive definite"


or the REML analysis is not converged with an error message

.. code:: sh

    "Log-likelihood not converged"


**Q1**: *How many samples are required for a GCTA-GREML analysis?*

A1: For unrelated individuals and common SNPs, you will need at least 3160 unrelated samples to get a SE down to 0.1 (see Visscher et al. 2014 PLoS Genet). For GREML analysis with multiple GRMs and/or GRM(s) computed from 1000G imputed data, a much larger sample size is required (see Yang et al. 2015 Nat Genet).

**Q2**: *Why do I need a small standard error (SE)?*

A2: The 95% confidence interval (CI) is approximately h2-SNP estimate +- 1.96 * SE. If the SE is too large, the 95% CI will cover the whole parameter space (from 0 to 1) so that you won't be able to make any meaningful inference from the estimate.


**2. How much memory do I need to run a GREML analysis?**

**1)** *Making a GRM*

This process involves a GRM, and a n x n matrix of the number of SNPs used for GRM calculation. 

Size of GRM in double :math:`precision = n * (n + 1) / 2 * 8` bytes

n x n matrix for the number of SNPs used to calculate GRM in single :math:`precision = n * (n + 1) /2 * 4` bytes

Thus the total memory is :math:`[n * (n + 1) / 2 * 12] / 10243 GB + 0.5GB`.

If the sample size n is huge, you can use --make-grm-part to reduce the memory usage (minimize and boost performance by parts divided). See :code:`--make-grm-part` for more details.

**2)** *REML analysis*

The REML process is a bit complicated. It involves a number of n x n matrices, e.g. GRM, variance-covariance V matrix, the projection P matrix and temporary matrices for V inverse calculation.

Total memory :math:`usage ~= (t + 4) * n * n * 8` bytes, where t is the number of genetic components (i.e. the number of GRMs) fitted in the model.

.. note::

    Note that these calculations haven't taken into account vectors and the other matrices of smaller size. Therefore, to submit a job to a computer cluster I would request 20% more memory than the predicted amount.

**3. How to calculate LRT in GREML?**

If there is only one genetic variance component (i.e. a single GRM) in your analysis, GCTA will calculate the LRT for the genetic variance automatically. The log likelihood for the full model (logL) and that for the reduced model (logL0) as well as the LRT and p-value will be reported in the \*.hsq file, where :math:`LRT = 2[logL - logL0]` which is distributed as a mixture of 0 and chi-squared (df = 1) with a probability of 0.5.

If you have multiple genetic variance components involved in your analysis (e.g. an analysis of genotype-environment (GE) interaction or a joint analysis of all chromosomes), by default, GCTA will only provide the LRT for first genetic variance component. In this case, you may need use the option :code:`--reml-lrt` to specify which component(s) you want to test. For example, for a GE interaction model, :math:`y = Xb + e + g + ge + e`, if you want to test the significance of the variance of GE interaction effects, you can add the option :math:`--reml-lrt 2` to your REML analysis:

.. code:: sh

    gcta64 --grm test --pheno test.phen --gxe test.gxe --reml --reml-lrt 2--out test


You can also calculate the LRT for multiple genetic variance components. For example, for a joint analysis of 22 chromosomes (22 genetic components in the model), you could test whether, for example, chromosomes 3 and 7 simultaneously by adding the option :code:`--reml-lrt 3 7` to the analysis:

You can also calculate the LRT for multiple genetic variance components. For example, for a joint analysis of 22 chromosomes (22 genetic components in the model), you could test whether, for example, chromosomes 3 and 7 simultaneously by adding the option :code:`--reml-lrt 3 7` to the analysis:

.. code:: sh

    gcta64 --mgrm grm_chrs.txt --pheno test.phen --reml --reml-lrt 3 7 --out test_chrs

The LRT for multiple components is distributed as a mixture of 0 and chi-squared (df = p) with a probability of 0.5, where p is the number of components to be tested.


**4. What does it mean if I get the following error messages?**

In MS Windows:

.. code::

    This application has requested the Runtime to terminate it in an unusual way.
    Please contact the application's support team for more information.

In Linux:

.. code:: sh
    
    terminate called after throwing an instance of 'std::bad_alloc'
      what():  St9bad_alloc
    Aborted

It means that the analysis requires more than 4 GB memory but the 32-bit version of GCTA only allows you to use a maximum of 4 GB memory. Solution: use the 64-bit version of GCTA on a 64-bit machine.


**5. Can I use GCTA in other species such as dogs and cattle?**

Yes, you can. You just need to specify the number of autosomes using the option --autosome-num when creating the GRM. For example:

.. code:: sh

    gcta64 --bfile test_dog --autosome-num 38 --autosome --make-grm --out test_dog

or

.. code:: sh

    gcta64 --bfile test_dog --autosome-num 38 --chr 1 --make-grm --out test_dog_c1
    gcta64 --bfile test_dog --autosome-num 38 --chr 2 --make-grm --out test_dog_c2
    ...
    gcta64 --bfile test_dog --autosome-num 38 --chr 38 --make-grm --out test_dog_c38

or

.. code:: sh

    gcta64 --bfile test_dog --autosome-num 38 --make-grm-xchr --out test_dog_xchr

**6. What does it mean if I get an estimate of V(G)/Vp of 0.9999?**

For a case-control study, V(G), V(e), Vp, V(G)/Vp are all on the observed scale. V(G)/Vp_L is the estimate of variance explained on the underlying liability scale under a threshold model. On the observed scale (0-1 disease status), the genetic variance can be greater Vp per definition, i.e. if the heritability on the underlying scale (h2L) is high and the disease prevalence is low, it is possible that the heritability on the observed scale (h2O) can be greater than 1. By default, GCTA does not allow any estimate of variance component to be negative. In this case, Ve is constrained at 10-6, so that the estimate of V(G)/Vp is constrained at 0.9999. You could specify the option --reml-no-constrain to allow V(G)/Vp to be greater than 1. However, you need to be cautious that any artefacts between cases and control will be estimated as 'genetic' variance, especially when cases and controls were genotyped separately (e.g. on different plate or at different labs). When using GCTA to analysis a case-control study, very stringent QC on SNPs are required. Please refer to Lee et al (2011 AJHG) for the QC steps and some other technical details of applying the method in case-control studies.

For a quantitative trait (which is relatively robust to the artefacts in SNP data as compared to a case-control study), it is likely that your sample size is small so that the estimate varies within a great range (i.e. large standard error). It may also suggest that the true parameter (i.e. variance explained by all SNPs) is relatively large.

**7. Can I use GCTA-GREML in family data?**

Yes, you can. GCTA-GREML does not assume that the individuals should be unrelated. The reason for excluding close-relatives in Yang et al. (Nat. Genet. 2010 and 2011) is because we do not want our estimates to be confounded with some possible shared environment effects and the effects of some possible causal variants that are not tagged by the SNPs but captured by pedigree information. If you are interested in the variance explained by a subset of SNPs in family data, you could fit the genetic relationship matrix (GRM) estimated from these SNPs along with a matrix of pedigree structure using the option --mgrm when running the REML analysis (--reml). Alternatively, we could fit the GRM of the subset of SNPs together with another GRM estimated from the SNPs in the rest of the genome.

If you don’t have SNP genotypes in the data and you are only interested in estimating pedigree-based heritability (see Yang et al. 2017 Nat Genet for definitions), we can compute a pedigree relatedness matrix from pedigree data using the script available here and run a REML analysis using the pedigree relatedness matrix as if it’s a GRM.

See GCTA-GREML in family data for an analysis of estimating SNP-based and pedigree-based h2 simultaneously in family data.


**8. Meta-analysis of GREML results from multiple cohorts**

If there are multiple cohorts and for some reason you are unable to pool all the individual-level genotype data together for a combined analysis, then it is OK to run a inverse-variance meta-analysis, i.e.

.. math::

    h^2meta = sum(h^2i / SE^2i) / sum(1 / SE^2i) with SE = (1 / sum(1 / SE^2i))^2

However, this is less powerful than a combined analysis because the meta-analysis does not utilise the contrasts between individuals across cohorts.


**9. Can I run a GREML analysis using a subset of SNPs selected by p-values from GWAS?**

If the SNPs are ascertained by p-value from GWAS analysis in the same sample, the GREML estimate of variance explained by this subset of SNPs will be inflated due to the winners' curse issue, i.e. the selection creates a positive correlation between true SNP effects and estimation errors.

If the SNPs are selected by p-values from association analysis in an independent sample, then it's OK. For example, in Wood et al. 2014 Nat Genet, we selected SNPs in a discovery set and performed GREML analysis of the selected SNPs in an independent validation set.


**10. Can I use the GRM to check for cryptic relatedness in my sample?**

Yes, you can. The expected value of Ajk = 

1\) 1 for MZ twins / duplicated samples

2\) 0.5 for 1st degree relatives (e.g. full-sibs or parent-offspring)

3\) 0.25 for 2nd degree relatives (e.g. grandparent-grandchild)

4\) 0.125 for 3rd degree relatives (e.g. cousins)

.. note::
    Note that these are the expected values. The realised GRM values come with sampling errors which is proportional to the number of markers used to compute the GRM. For distant relatives (e.g. cousins 2 times removed), we might not have enough power (or precision) distinguish them from unrelated pairs. See Supplementary Note #2 of Yang et al. (2010 Nature Genetics) for more details.

There are two ways of reading the GRM in R.

* See the sample code for reading the binary GRM file.

* Using --make-grm-gz option to convert the binary format to compressed text format.


**11. Can I run a GBLUP prediction analysis with GCTA?**

1) Creating a GRM using SNP data

.. code:: sh

    gcta64  --bfile test  --make-grm test  --out test

2) REML analysis with the --reml-pred-rand option to output the BLUP solutions of the individuals (i.e. estimate of total genetic value of each individual)

.. code:: sh

    gcta64 --reml --grm test --pheno test.phen --reml-pred-rand --out test


From the analysis above, you will have a output file test.indi.blp. There is no header line. Columns are family ID, individual ID, an intermediate variable, the total genetic value, another intermediate variable and the residual. If there are multiple GRMs fitted in the REML analysis, each GRM will insert additional two columns, i.e. an intermediate variable and a total genetic value, in front of the last two columns.

::

    01       0101    -0.012    -0.014   -0.010    -0.035
    02       0203    0.021     0.031    -0.027    -0.031
    03       0305    0.097     0.102    -0.026    -0.041

For a mixed linear model y = g + e, the BLUP estimates of genetic values (ug) and residuals (ue) are calculated using the two equations below (Lynch and Walsh 1996, page 749)

.. math::

    g_{hat} = V_gAV^{-1}y\ and\ e_{hat} = V_eV^{-1}y

where Vg is the genetic variance, Ve is the residual variance, A is the GRM, and y is the phenotype vector.

3\) BLUP solutions for the SNP effects

.. code:: sh

    gcta64 --bfile test --blup-snp test.indi.blp --out test

The result will be saved in a file test.snp.blp. Columns are SNP ID, reference allele and BLUP of SNP effect. If there are multiple GRMs, each GRM will add an additional column to the file. You can alway ignore the last column.


::
    
    rs103645   A     0.00312    0.00451
    rs175292   G     -0.00021   0.00139


4\) You may then use PLINK --score option using the test.snp.blp as input to predict the polygenic profiles of new samples.


**12. Can I run a bivariate GCTA-GREML of two independent samples?**

Bivariate GCTA-GREML of two independent samples

Here is an example of performing a bivariate GCTA-GREML analysis for two traits measured in two independent samples.

1\) Creating a GRM for all the individuals combined (from the two samples)

2\) Creating a phenotype file of two traits for all the samples. Assuming 100 individuals in sample #1 and 100 individuals in sample #2, here is an example of the phenotype file ("NA" represents missing data)

::

    FID    IID    trait1   trait2
    1      1      0.1      NA
    2      2      0.2      NA
    3      3      0.1      NA
    ...
    100    100    0.5      NA
    101    101    NA       2.1
    102    102    NA       3.1
    103    103    NA       2.2
    ...
    200    200    NA       2.1

.. note::
    
    this analysis also applies to a single trait measured in two samples. Then the analysis is to estimate genetic correlation between two samples for the same trait.

**13. How can I estimate the fixed effects from GCTA-GREML?**

For an analysis without a covariate, the GREML model can be written as:

.. math::

    y = mu + g + e

where mu is the mean term (fixed effect), g is the genetic value (random effect) and e is the residual.

1\) Categorical covariate (e.g. sex and cohort): :code:`--covar` option

If the covariate is a categorical covariate, there will be t - 1 variables (where t is the number of categories, e.g. t = 2 for sex) because otherwise the XTV-1X will not be invertible (X is design matrix for the fixed effects and V is the covariance-covariance matrix). Therefore, the model can be written as

.. math::
    
    y = mu + x_{c(2)}*b_{c(2)} + x_{c(3)}*b_{c(3)} + … + x_{c(t)}*b_{c(t)} + g + e

where x is coded as 1 or 0 (representing the presence or absence of a category), bc(i) is interpreted as difference in mean phenotype in category i from the category 1. Note that the order of the categories are determined by their order of appearance in the data.

2\) Quantitative covariate (e.g. age): :code:`--qcovar` option


The covariate is fitted as a continuous variable, then the model is :math:`y = mu + x_{q(1)}*b_{q(1)}` + g + e where the interpretation of bq(1) is similar as that from a linear regression.

3\) If we have a categorical covariate and two quantitative covariates, the model is

.. math::
    
    y = mu + x_{c(2)}*b_{c(2)} + x_{c(3)}*b_{c(3)} + … + x_{c(t)}*b_{c(t)} + x_{q(1)}*b_{q(1)} + x_{q(2)}*b_{q(2)} + g + e

Of course, we could also fit multiple quantitative covariates and multiple categorical covariates.

.. note::
    
    These fixed effects can be estimated using the :code:`--reml-est-fix` option in a REML analysis. The estimates are shown in the log output following the order in the model above, i.e. the effect of each quantitative covariate followed by the effect each of category of the categorical covariates.

**14. Why do I get a negative estimate of SNP-heritability?**

Heritability (:math:`h^2`) is per definition non-negative. However, the estimate of :math:`h^2` is supposed to be following a normal distribution with mean :math:`h^2` and variance :math:`SE^2` where SE is the standard error of the estimate of :math:`h^2`. Therefore, to get an unbiased estimate of :math:`h^2`, we should allow the estimate to be negative (:code:`--reml-no-constrain` option in GCTA-GREML analysis).

In practice, there are a least two scenarios when we would see negative estimate of :math:`h^2`.

* Small sample size. If the sample size is small, the sampling variance (:math:`SE^2`) will be large. In this case, the estimate of :math:`h^2` will fluctuate a lot and therefore has a certain chance to jump out of the parameter space (between 0 and 1).

* The true :math:`h^2` parameter is small. If :math:`h^2` is very small, then even if the sample size is large, we will still have a certain probability to see negative estimate.

In the Yang et al. (2013 PLoS Genet) and Zhu et al. (2015 AJHG) papers, to get an unbiased estimate of the mean estimate of :math:`h^2`, we did not constrain the estimate to 0.

**15. Error: variance-covaraince matrix V is not positive definite**

The GREML method uses REML for variance estimation (please see Yang et al. 2010 AJHG for details), which requires the inverse of the variance-covariance matrix V. If V is not positive definite, the inverse of V does not exist. We therefore could not estimate the variance component. This usually happens when one (or more) of the variance components are negative or constrained at zero. It might also indicate there is something wrong with the GRM or the data which you might need to check carefully.

Unfortunately, there has not been an ultimate solution. Tricks such as adding a small number of to the diagonal elements of V also do not guarantee the modified V being invertible. In some cases, you might be able to get around the problem by using alternative REML algorithms e.g. the Fisher scoring approach (--reml-alg 1).

We have implemented the "bending" approach (Hayes and Hill 1981 Biometrics) in GCTA to invert V if V is not positive definite (you could add the --reml-bendV option to a REML or MLMA analysis to activate this approach). The "bending" approach guarantees to get an approximate of V-1 but it does not guarantee the REML analysis being converged.

.. note::

    Note that the :code:`--reml-bendV` option only provides an approximate inverse of V and has not been tested extensively. The results from analyses using this option might not be reliable.

**16. GREML p-value = 0?**

This is a precision issue. It means that the p-value is extremely small. You can calculate a more precise p-value in R.

1) :code:`p-value = 0.5 * pchisq(LRT, df=1, lower.tail=FALSE) # one-tailed test`, e.g. :math:`h^2_g` is constrained to be positive in a GREML analysis.

2) :code:`p-value = pchisq(LRT, df=1, lower.tail=FALSE) # two-tailed test` (recommended to test whether :math:`r_g = 0` in a bivariate GREML analysis or to test if :math:`h^2_g = 0` in a unconstrained GREML analysis).

No LRT reported in \*.hsq output file?

.. math::

    LRT ~= (estimate / SE)^2
