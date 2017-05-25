# monocle-wrapper-for-Galaxy
Monocle wrapper for Galaxy
https://bioconductor.org/packages/release/bioc/html/monocle.html

Creating the CellDataSet class

From Monocle Vignette and Documentation:
The monocle package takes a matrix of gene expression values as calculated by Cuffinks [2] or another gene expression estimation program. Monocle can work with relative expression values (e.g. FPKM or TPM units) or absolute transcript counts (e.g. from UMI experiments). Monocle also works "out-of-the-box" with the transcript count matrices produced by CellRanger, the software pipeline for analyzing experiments from the 10X Genomics Chromium instrument. Monocle also works well with data from other RNA-Seq work flows such as sci-RNA-Seq and instruments like the Biorad ddSEQ.<br>
Although Monocle can be used with raw read counts, these are not directly proportional to expression values unless you normalize them by length, so some Monocle functions could produce nonsense results. If you don't have UMI counts, We recommend you load up FPKM or TPM values instead of raw read counts.

*monocle* holds single cell expression data in objects of the *CellDataSet* class.  The class is derived from the Bioconductor *ExpressionSet* class, which provides a common interface familiar to those who have analyzed microarray experiments with Bioconductor.  The class requires three input files:
<br>
1. exprs, a numeric matrix of expression values, where rows are genes, and columns are cells

2. phenoData, an *AnnotatedDataFrame* object, where rows are cells, and columns are vell attributes (such as cell type, culture condition, day captured, etc.)

3. featureData, an *AnnotatedDataFrame* object, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.

Function(s)
------
### newCellDataSet
*Creates a new CellDataSet object*
<dl>
	<strong>Arguments</strong>
	<dt>cellData</dt>
	<dd>expression data matrix for an experiment</dd>
	<dt>phenoData</dt>
	<dd>data frame containing attributes of individual cells</dd>
	<dt>featureData</dt>
	<dd>data frame containing attributes of features (e.g. genes)</dd>
	<dt>lowerDetectionLimit</dt>
	<dd>the minimum expression level that consistitutes true expression</dd>
	<dt>expressionFamily</dt>
	<dd>the VGAM family function to be used for expression response variables</dd>
</dl>

Notes - WIP:
- show mRNA totals