# VcfEigen01

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotator for the data of https://xioniti01.u.hpc.mssm.edu/v1.1/ : Eigen makes use of a variety of functional annotations in both coding and noncoding regions (such as made available by the ENCODE and Roadmap Epigenomics projects), and combines them into one single measure of functional importance.


## Usage

```
Usage: vcfeigen [options] Files
  Options:
  * -D, --directory
      Eigen directory containing the tabix files *.tab.gz
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -p, --prefix, --tabixFilePrefix
      prefix of the files in the eigen directory
      Default: Eigen_hg19_
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * variant


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfeigen
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfeigen/VcfEigen01.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfeigen/VcfEigen01.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfeigen** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### About Eigen

"Eigen is a spectral approach to the functional annotation of genetic variants in coding and noncoding regions. Eigen makes use of a variety of functional annotations in both coding and noncoding regions (such as made available by the ENCODE and Roadmap Epigenomics projects), and combines them into one single measure of functional importance. Eigen is an unsupervised approach, and, unlike most existing methods, is not based on any labelled training data. Eigen produces estimates of predictive accuracy for each functional annotation score, and subsequently uses these estimates of accuracy to derive the aggregate functional score for variants of interest as a weighted linear combination of individual annotations. We show that the resulting meta-score has good discriminatory ability using disease associated and putatively benign variants from published studies (for both Mendelian and complex diseases). The Eigen score is particularly useful in prioritizing likely causal variants in a region of interest when it is combined with population-level genetic data in the framework of a hierarchical model. Furthermore, an important advantage of the Eigen score is that it can be easily adapted to a specific tissue or cell type. More information about the Eigen score can be found in the accompanying manuscript: [A spectral approach integrating functional genomic annotations for coding and noncoding variants (Iuliana Ionita-Laza, Kenneth McCallum, Bin Xu, Joseph Buxbaum).A spectral approach integrating functional genomic annotations for coding and noncoding variants (Iuliana Ionita-Laza, Kenneth McCallum, Bin Xu, Joseph Buxbaum).](http://www.columbia.edu/~ii2135/Eigen_11_24.pdf) 

### Example

```

$ gunzip -c input.vcf.gz | cut -f 1-8 | java -jar dist/vcfeigen01.jar -D eigen_data_dir  | grep -i eigen

##INFO=<ID=EIGEN_CODING_Consequence,Number=A,Type=String,Description="Consequence">
##INFO=<ID=EIGEN_CODING_Eigen_PC_phred,Number=A,Type=Float,Description="Eigen-PC-phred">
##INFO=<ID=EIGEN_CODING_Eigen_PC_raw,Number=A,Type=Float,Description="Eigen-PC-raw">
##INFO=<ID=EIGEN_CODING_Eigen_phred,Number=A,Type=Float,Description="Eigen-phred">
##INFO=<ID=EIGEN_CODING_Eigen_raw,Number=A,Type=Float,Description="Eigen-raw">
##INFO=<ID=EIGEN_CODING_GERP_NR,Number=A,Type=Float,Description="GERP_NR">
##INFO=<ID=EIGEN_CODING_GERP_RS,Number=A,Type=Float,Description="GERP_RS">
##INFO=<ID=EIGEN_CODING_MA,Number=A,Type=Float,Description="MA">
##INFO=<ID=EIGEN_CODING_PhastPla,Number=A,Type=Float,Description="PhastPla">
##INFO=<ID=EIGEN_CODING_PhastPri,Number=A,Type=Float,Description="PhastPri">
##INFO=<ID=EIGEN_CODING_PhastVer,Number=A,Type=Float,Description="PhastVer">
##INFO=<ID=EIGEN_CODING_PhyloPla,Number=A,Type=Float,Description="PhyloPla">
##INFO=<ID=EIGEN_CODING_PhyloPri,Number=A,Type=Float,Description="PhyloPri">
##INFO=<ID=EIGEN_CODING_PhyloVer,Number=A,Type=Float,Description="PhyloVer">
##INFO=<ID=EIGEN_CODING_PolyPhenDIV,Number=A,Type=Float,Description="PolyPhenDIV">
##INFO=<ID=EIGEN_CODING_PolyPhenVar,Number=A,Type=Float,Description="PolyPhenVar">
##INFO=<ID=EIGEN_CODING_SIFT,Number=A,Type=Float,Description="SIFT">
##INFO=<ID=EIGEN_NC_DnasePval,Number=A,Type=Float,Description="DnasePval">
##INFO=<ID=EIGEN_NC_DnaseSig,Number=A,Type=Float,Description="DnaseSig">
##INFO=<ID=EIGEN_NC_Eigen_PC_phred,Number=A,Type=Float,Description="Eigen-PC-phred">
##INFO=<ID=EIGEN_NC_Eigen_PC_raw,Number=A,Type=Float,Description="Eigen-PC-raw">
##INFO=<ID=EIGEN_NC_Eigen_phred,Number=A,Type=Float,Description="Eigen-phred">
##INFO=<ID=EIGEN_NC_Eigen_raw,Number=A,Type=Float,Description="Eigen-raw">
##INFO=<ID=EIGEN_NC_FairePval,Number=A,Type=Float,Description="FairePval">
##INFO=<ID=EIGEN_NC_FaireSig,Number=A,Type=Float,Description="FaireSig">
##INFO=<ID=EIGEN_NC_GERP_NR,Number=A,Type=Float,Description="GERP_NR">
##INFO=<ID=EIGEN_NC_GERP_RS,Number=A,Type=Float,Description="GERP_RS">
##INFO=<ID=EIGEN_NC_H3K27ac,Number=A,Type=Float,Description="H3K27ac">
##INFO=<ID=EIGEN_NC_H3K4Me1,Number=A,Type=Float,Description="H3K4Me1">
##INFO=<ID=EIGEN_NC_H3K4Me3,Number=A,Type=Float,Description="H3K4Me3">
##INFO=<ID=EIGEN_NC_OCPval,Number=A,Type=Float,Description="OCPval">
##INFO=<ID=EIGEN_NC_PhastPla,Number=A,Type=Float,Description="PhastPla">
##INFO=<ID=EIGEN_NC_PhastPri,Number=A,Type=Float,Description="PhastPri">
##INFO=<ID=EIGEN_NC_PhastVer,Number=A,Type=Float,Description="PhastVer">
##INFO=<ID=EIGEN_NC_PhyloPla,Number=A,Type=Float,Description="PhyloPla">
##INFO=<ID=EIGEN_NC_PhyloPri,Number=A,Type=Float,Description="PhyloPri">
##INFO=<ID=EIGEN_NC_PhyloVer,Number=A,Type=Float,Description="PhyloVer">
##INFO=<ID=EIGEN_NC_PolIIPval,Number=A,Type=Float,Description="PolIIPval">
##INFO=<ID=EIGEN_NC_PolIISig,Number=A,Type=Float,Description="PolIISig">
##INFO=<ID=EIGEN_NC_TFBS_max,Number=A,Type=Float,Description="TFBS_max">
##INFO=<ID=EIGEN_NC_TFBS_num,Number=A,Type=Float,Description="TFBS_num">
##INFO=<ID=EIGEN_NC_TFBS_sum,Number=A,Type=Float,Description="TFBS_sum">
##INFO=<ID=EIGEN_NC_cmycPval,Number=A,Type=Float,Description="cmycPval">
##INFO=<ID=EIGEN_NC_cmycSig,Number=A,Type=Float,Description="cmycSig">
##INFO=<ID=EIGEN_NC_ctcfPval,Number=A,Type=Float,Description="ctcfPval">
##INFO=<ID=EIGEN_NC_ctcfSig,Number=A,Type=Float,Description="ctcfSig">
(...)
14	741	.	A	C	1	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0;EIGEN_NC_Eigen_PC_phred=7.69892;EIGEN_NC_Eigen_PC_raw=-0.05363416;EIGEN_NC_Eigen_phred=1.64761;EIGEN_NC_Eigen_raw=-0.27013662;EIGEN_NC_FairePval=0.0;EIGEN_NC_FaireSig=0.0;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=62.64;EIGEN_NC_H3K4Me1=9.0;EIGEN_NC_H3K4Me3=13.04;EIGEN_NC_OCPval=0.0;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.0;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0;SNVHPOL=4;SNVSB=0.0
14	1142	.	T	A	19	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=13.9238;EIGEN_NC_Eigen_PC_raw=0.7665708;EIGEN_NC_Eigen_phred=3.47555;EIGEN_NC_Eigen_raw=-0.10196207;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=56.68;EIGEN_NC_H3K4Me1=7.48;EIGEN_NC_H3K4Me3=13.0;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=4;SNVSB=0.0
14	1195	.	T	A	2	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=13.9238;EIGEN_NC_Eigen_PC_raw=0.7665708;EIGEN_NC_Eigen_phred=3.47555;EIGEN_NC_Eigen_raw=-0.10196207;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=56.68;EIGEN_NC_H3K4Me1=7.48;EIGEN_NC_H3K4Me3=13.0;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=4;SNVSB=0.0
14	1437	.	T	A	1	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=13.9751;EIGEN_NC_Eigen_PC_raw=0.7803033;EIGEN_NC_Eigen_phred=3.53009;EIGEN_NC_Eigen_raw=-0.09754433;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=61.56;EIGEN_NC_H3K4Me1=7.64;EIGEN_NC_H3K4Me3=14.64;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=4;SNVSB=0.0
14	1659	.	G	C	8	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=14.1281;EIGEN_NC_Eigen_PC_raw=0.82222354;EIGEN_NC_Eigen_phred=3.69067;EIGEN_NC_Eigen_raw=-0.08474564;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=78.24;EIGEN_NC_H3K4Me1=7.56;EIGEN_NC_H3K4Me3=18.84;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=3;SNVSB=0.0
14	1910	.	T	C	58	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=13.8695;EIGEN_NC_Eigen_PC_raw=0.75230044;EIGEN_NC_Eigen_phred=3.41565;EIGEN_NC_Eigen_raw=-0.106857374;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=52.16;EIGEN_NC_H3K4Me1=5.84;EIGEN_NC_H3K4Me3=12.16;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=4;SNVSB=-10.8


```

