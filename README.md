# TRACS
* Prerequisite : R (>=3.2.0, dtwclust package), Python 2.* (Scipy, numpy, matplotlib, sklearn, rpy2).
* TRACS takes six (+ three optional) arguments
* Example:
```console
python TRACS.py -f expr.txt -tn 5 -rn 3 -tp 0,10,20,40,60 -o results -og hsa
```

## Arguments
1. Gene expression file (-f, --file)
   * Gene expression file should be tab-delimited.
   * If there are T time points for each biological replicate, columns should be arranges as <br>[T columns from first replicate] -> [T columns from second replicate] -> ...
   * Example:
```
Gene  R1-T1  R1-T2 R1-T3 ... R2-T1 R2-T2 R2-T3 ...
Gene1 8 11  8 ... 19  10  11  ...
... 19  21  ... 11  10  11  ...
Gene2 14  20  ... 10  8 20  ...
```

2. The number of time points (-tn, --timenums) and replicates (-rn, --repnums)
3. A list of time points (-tp, --timepoints)
   * A list of time points is the format of t1,t2,t3,...,tn
   * Example: ``` -tn 0,5,10,15,20,25  ```

4. Output directory (-o, --outdir)
   * All the output files will be generated in the output directory

5. Organism of interest (-og, --organism)
   * KEGG organism code (https://www.kegg.jp/kegg/catalog/org_list.html)
   * Example: when your data is human gene expression data ``` -og hsa ```
 
6. Method for clustering (-m, --method) (optional)
   * Select among KM (K-means clustering), AC (Agglomerative clustering), KS (K-Shape)
   * Default m=KM

7. Start (-ks, --kstart) end (-ke, --kend) of the number of clusters, K (optional)
   * TRACS will search for the optimal K in the range of [ks, ke].
   * Default ks=1, ke=10

