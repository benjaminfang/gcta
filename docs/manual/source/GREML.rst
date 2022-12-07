GREML
===========

Tutorial
---------------

If you have used PLINK before, you will find it easy to use GCTA. In this tutorial, all the options used are not detailed. Please refer to the documentation of GCTA for details of the options and formats of the input or output files.

**GCTA-GRM: calculating the genetic relationship matrix (GRM) from all the autosomal SNPs**

Suppose you have a GWAS data set in PLINK binary PED format, e.g. test.bed, test.bim and test.fam. You can type this command to calculate the genetic relationships between pairwise individuals from all the autosomal SNPs

.. code:: sh
    
    gcta64 --bfile test --autosome --maf 0.01 --make-grm --out test --thread-num 10

The genetic relationship matrix will be saved in the files test.grm.bin, test.grm.N.bin and test.grm.id .

For datasets with an extremely large number of SNPs and large sample size (e.g. 1000G imputed data, you can use the following commands:

.. code:: sh

    gcta64 --bfile test --chr 1 --maf 0.01 --make-grm --out test_chr1 --thread-num 10
    gcta64 --bfile test --chr 2 --maf 0.01 --make-grm --out test_chr2 --thread-num 10
    ...
    gcta64 --bfile test --chr 22 --maf 0.01 --make-grm --out test_chr22 --thread-num 10

which calculate the GRM for each autosome and then merge the 22 GRMs by the following command:

.. code:: sh
    
    gcta64 --mgrm grm_chrs.txt --make-grm --out test

You can use this command to remove cryptic relatedness:

.. code:: sh

    gcta64 --grm test --grm-cutoff 0.025 --make-grm --out test_rm025

which creates a new GRM of "unrelated" individuals. Please be aware that the cutoff value 0.025 is quite arbitrary.


**GCTA-GREML analysis: estimating the variance explained by the SNPs**

.. code:: sh

    gcta64 --grm test --pheno test.phen --reml --out test --thread-num 10

The results will be saved in the file test.hsq. 

You can also include the first 4 or 10 eigenvectos from principal component analysis (PCA) as covariates by the command

.. code:: sh

    gcta64 --grm test --pheno test.phen --reml --qcovar test_10PCs.txt --out test --thread-num 10

You can also estimate the variance explained by the SNPs on each chromosome by fitting one chromosome at a time:

.. code:: sh

    gcta64 --grm test_chr1 --pheno test.phen --reml --out test_chr1 --thread-num 10
    gcta64 --grm test_chr2 --pheno test.phen --reml --out test_chr2 --thread-num 10
    ......
    gcta64 --grm test_chr22 --pheno test.phen --reml --out test_chr22 --thread-num 10

or fitting all the 22 autosomes simultaneously by

.. code:: sh

    gcta64 --mgrm grm_chrs.txt --pheno test.phen --reml --out test_all_chrs --thread-num 10

You are also allowed to include the first 4 or 10 eigenvectors from PCA as covariates in any of these analyses.

**GCTA-GREML analysis for a case-control study**

For a case-control study, the phenotypic values of cases and controls should be specified as 1 and 0, respectively. Suppose you have prepared a phenotype file test_cc.phen. You can type the following command to estimate the variance explained by all the autosomal SNPs on the observed 0-1 scale and transform the estimate to that on the underlying liability scale (assuming the disease prevalence is 0.01 in this example)

.. code:: sh

    gcta64 --grm test --pheno test_cc.phen --reml --prevalence 0.01 --out test --thread-num 10


Making a GRM
---------------

**GCTA-GRM: estimating genetic relatedness from SNPs**

``--make-grm`` or ``--make-grm-bin``
    Estimate the genetic relationship matrix (GRM) between pairs of individuals from a set of SNPs and save the lower triangle elements of the GRM to binary files, e.g. test.grm.bin, test.grm.N.bin, test.grm.id.

    Output file
        test.grm.bin (it is a binary file which contains the lower triangle elements of the GRM).
        test.grm.N.bin (it is a binary file which contains the number of SNPs used to calculate the GRM).
        test.grm.id (no header line; columns are family ID and individual ID, see above).
        You can not open test.grm.bin or test.grm.N.bin by a text editor but you can use the following R script to read them in R)

    .. code:: R
        
        # R script to read the GRM binary file
        ReadGRMBin=function(prefix, AllN=F, size=4){
          sum_i=function(i){
            return(sum(1:i))
          }
          BinFileName=paste(prefix,".grm.bin",sep="")
          NFileName=paste(prefix,".grm.N.bin",sep="")
          IDFileName=paste(prefix,".grm.id",sep="")
          id = read.table(IDFileName)
          n=dim(id)[1]
          BinFile=file(BinFileName, "rb");
          grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
          NFile=file(NFileName, "rb");
          if(AllN==T){
            N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
          }
          else N=readBin(NFile, n=1, what=numeric(0), size=size)
          i=sapply(1:n, sum_i)
          return(list(diag=grm[i], off=grm[-i], id=id, N=N))
        }
            

    .. note::
        
        --make-grm has been rewritten with orders of magnitude improvement in speed and memory usage. Currently, It can only used in combination with a limited number of other flags, i.e., --keep, --remove, --chr, --autosome-num, --autosome, --extract, --exclude, --maf, --max-maf, --thread-num, --update-ref-allele, --update-sex, --update-freq. You can use --make-grm-part to reduce the memory usage further.

    Make GRM function can combine with --mbfile to calculate GRMs in multiple PLINK files without merge them together. 

--mbfile chrs.txt
    If the genotype data is very large, the data is often saved in separate PLINK files (e.g. one for each chromosome). Use --mbfile to specify multiple PLINK files. The input is a text file with each row representing a PLINK binary file (without file name suffix).

    Input file format: ::
        
        data_chr1
        data_chr2
        …

    .. note::
        
        All these files shall have same sample size and order, the program will prompt an error if not. 

--make-grm-part m i
    Partition the GRM into m parts (by row), and compute the i-th part in the current run. 

    .. note::
        
        This option is designed to compute the GRM in a very large sample (e.g. the UK Biobank data). The memory usage of each run is the total memory required divided by m. Thus partitioning a large number of parts can reduce the memory usage significantly. The total memory required is approximately [n * (n + 1) / 2 * 12] / 10243 GB + 0.5GB, where n is the sample size. As some computer clusters limit the virtual memory, allocating 1 to 2GB more memory to each job will be safer. In our computation of the GRM in the UKB data, we partitioned the whole data set (n = 456,426) into 250 parts and allocated 6700MB memory to each job.


    Example:

    .. code:: sh
        
        # Partition the GRM into 3 parts
        gcta64 --bfile test --make-grm-part 3 1 --thread-num 5 --out test
        gcta64 --bfile test --make-grm-part 3 2 --thread-num 5 --out test
        gcta64 --bfile test --make-grm-part 3 3 --thread-num 5 --out test
        # Merge all the parts together (Linux, Mac)
        cat test.part_3_*.grm.id > test.grm.id
        cat test.part_3_*.grm.bin > test.grm.bin
        cat test.part_3_*.grm.N.bin > test.grm.N.bin
        # Windows alternative
        copy /b test.part_3_*.grm.id test.grm.id
        copy /b test.part_3_*.grm.bin test.grm.bin
        copy /b test.part_3_*.grm.N.bin test.grm.N.bin
        


--make-grm-alg 0
    The default value is 0, and the GRM is calculated using the equation sum{[(xij - 2pi)*(xik - 2pi)] / [2pi(1-pi)]} as described in Yang et al. 2010 Nat Genet. If the value = 1, the GRM will be calculated using the equation sum[(xij - 2pi)(xik - 2pi)] / sum[2pi(1-pi)]*.


--make-grm-gz
    Estimate the GRM, save the lower triangle elements to a compressed text file (e.g. test.grm.gz) and save the IDs in a plain text file (e.g. test.grm.id).

    Output file format.  test.grm.gz (no header line; columns are indices of pairs of individuals (row numbers of the test.grm.id), number of non-missing SNPs and the estimate of genetic relatedness) ::
        
        1    1    1000    1.0021
        2    1    998     0.0231
        2    2    999     0.9998
        3    1    1000    -0.0031

    test.grm.id (no header line; columns are family ID and individual ID) ::

        011      0101
        012      0102
        013      0103
        ...

--make-grm-xchr
    Estimate the GRM from SNPs on the X-chromosome. The GRM will be saved in the same binary format as above (\*.grm.bin, \*.grm.N.bin and \*.grm.id). Due to the speciality of the GRM for the X-chromosome, it is not recommended to manipulate the matrix by --grm-cutoff or --grm-adj, or merge it with the GRMs for autosomes (see below for the options of manipulating the GRM).

    .. note::
        
        This flag has been re-implemented in GCTA 1.91.4, it has same performance and memory consumption as --make-grm.

    .. note::
        
        The function treats X chr as non-pseudoautosomal region (nPAR) with genotype coding for male as 0, 2. For pseudoautosomal region (PAR), we can alter the chromosome number in bim file to autosome and use --make-grm to run. Don't put nPAR and PAR together as X chr, GCTA will give weird results.


--make-grm-xchr-part m i
    Partition the GRM of X chromosome into m parts (by row), and compute the i-th part in the current run.
    
    See the document of --make-grm-part

--make-grm-xchr-gz
    Same as --make-grm-xchr but the GRM will be in compressed text files (see --make-grm-gz for the format of the output files).

--make-grm-inbred or --make-grm-inbred-gz
    Make a GRM for an inbred population such as inbred mice or inbred crops.

--ibc
    Estimate the inbreeding coefficient from the SNPs by 3 different methods.

    Output file format.  test.ibc (one header line; columns are family ID, individual ID, number of nonmissing SNPs, estimator 1, estimator 2 and estimator 3) ::

        FID      IID        NOMISS       Fhat1       Fhat2         Fhat3
        011      0101       999          0.00210     0.00198       0.00229
        012      0102       1000         -0.0033     -0.0029       -0.0031
        013      0103       988          0.00120     0.00118       0.00134

    See Yang et al. 2011 AJHG for the definitions of Fhat1, Fhat2 and Fhat3.

    Examples:

    .. code:: sh

        # Estimate the GRM from all the autosomal SNPs
        gcta64  --bfile test  --autosome  --make-grm  --out test

        # Estimate the GRM from the SNPs on the X-chromosome
        gcta64  --bfile test  --make-grm-xchr  --out test_xchr

        # Estimate the GRM from the SNPs on chromosome 1 with MAF from 0.1 to 0.4
        gcta64  --bfile test  --chr 1  --maf 0.1  --max-maf 0.4  --make-grm  --out test

        # Estimate the GRM using a subset of individuals and a subset of autosomal SNPs with MAF < 0.01
        gcta64  --bfile test  --keep test.indi.list  --extract test.snp.list  --autosome  --maf 0.01 --make-grm  --out test

        # Estimate the GRM from the imputed dosage scores for the SNPs with MAF > 0.01 and imputation R2 > 0.3
        gcta64  --dosage-mach  test.mldose.gz  test.mlinfo.gz  --imput-rsq  0.3  --maf 0.01  --make-grm --out test

        # Estimate the GRM from the imputed dosage scores for a subset of individuals and a subset of SNPs
        gcta64  --dosage-mach  test.mldose.gz  test.mlinfo.gz  --keep test.indi.list  --extract test.snp.list  --make-grm --out test

        # Estimate the inbreeding coefficient from all the autosomal SNPs
        gcta64  --bfile test  --autosome  --ibc  --out test

        # Calculate the GRM using the alternative method
        gcta64  --bfile test  --autosome --make-grm  --make-grm-alg 1  --out test_alg1


**<span>Citations</span>**
    *Method for estimating the GRM*:
        Yang et al. (2010) Common SNPs explain a large proportion of the heritability for human height. Nat Genet. 42(7): 565-9. [PubMed ID: 20562875]

    *Method for estimating the inbreeding coefficients and GCTA software*:
        Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]
