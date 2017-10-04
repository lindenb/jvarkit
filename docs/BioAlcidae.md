# BioAlcidae

javascript version of awk for bioinformatics


## Usage

```
Usage: bioalcidae [options] Files
  Options:
    -e, --expression
      Javascript expression
    -F, --format
      force format: one of VCF BAM SAM FASTQ FASTA BLAST DBSNP. BLAST is BLAST 
      XML version 1. DBSNP is XML output of NCBI dbsnp. INSDSEQ is XML output 
      of NCBI EFetch rettype=gbc.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -J, --json
      Optional. Reads a JSON File using google gson 
      (https://google-gson.googlecode.com/svn/trunk/gson/docs/javadocs/index.html 
      ) and injects it as 'userData' in the javascript context.
    -o, --output
      Output file. Optional . Default: stdout
    -f, --scriptfile
      Javascript file
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * vcf
 * javascript
 * js
 * nashorn



## See also in Biostars

 * [https://www.biostars.org/p/257346](https://www.biostars.org/p/257346)
 * [https://www.biostars.org/p/183197](https://www.biostars.org/p/183197)
 * [https://www.biostars.org/p/185162](https://www.biostars.org/p/185162)
 * [https://www.biostars.org/p/153060](https://www.biostars.org/p/153060)
 * [https://www.biostars.org/p/152016](https://www.biostars.org/p/152016)
 * [https://www.biostars.org/p/152720](https://www.biostars.org/p/152720)
 * [https://www.biostars.org/p/152820](https://www.biostars.org/p/152820)
 * [https://www.biostars.org/p/218444](https://www.biostars.org/p/218444)
 * [https://www.biostars.org/p/224402](https://www.biostars.org/p/224402)
 * [https://www.biostars.org/p/241751](https://www.biostars.org/p/241751)
 * [https://www.biostars.org/p/240452](https://www.biostars.org/p/240452)
 * [https://www.biostars.org/p/248385](https://www.biostars.org/p/248385)
 * [https://www.biostars.org/p/186610](https://www.biostars.org/p/186610)
 * [https://www.biostars.org/p/242127](https://www.biostars.org/p/242127)
 * [https://www.biostars.org/p/167389](https://www.biostars.org/p/167389)
 * [https://www.biostars.org/p/187494](https://www.biostars.org/p/187494)
 * [https://www.biostars.org/p/183197](https://www.biostars.org/p/183197)
 * [https://www.biostars.org/p/152820](https://www.biostars.org/p/152820)
 * [https://www.biostars.org/p/178004](https://www.biostars.org/p/178004)
 * [https://www.biostars.org/p/156250](https://www.biostars.org/p/156250)
 * [https://www.biostars.org/p/202400](https://www.biostars.org/p/202400)
 * [https://www.biostars.org/p/183982](https://www.biostars.org/p/183982)
 * [https://www.biostars.org/p/173201](https://www.biostars.org/p/173201)


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make bioalcidae
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bioalcidae/BioAlcidae.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bioalcidae/BioAlcidae.java)


<details>
<summary>Git History</summary>

```
Fri Aug 4 16:40:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/57f08e720a97f952bab81961431d83accdefeae3
Wed Jun 14 17:01:36 2017 +0200 ; fast genotype gvcf ; https://github.com/lindenb/jvarkit/commit/d77e93940ad9a7f8144527332067b663b55a10f6
Thu May 18 10:38:01 2017 +0200 ; vcfilterjs ; https://github.com/lindenb/jvarkit/commit/437e95c1e567ec25812c12608f3df409f24c2fa1
Wed May 17 18:02:11 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/7d5989fdf77a7a4a616a29876997809683224316
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Fri May 12 18:07:46 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/ca96bce803826964a65de33455e5231ffa6ea9bd
Sun May 7 13:21:47 2017 +0200 ; rm xml ; https://github.com/lindenb/jvarkit/commit/f37088a9651fa301c024ff5566534162bed8753d
Thu Apr 20 17:17:22 2017 +0200 ; continue transition jcommander ; https://github.com/lindenb/jvarkit/commit/fcf5def101925bea9ddd001d8260cf65aa52d6a0
Fri Mar 17 17:06:05 2017 +0100 ; added SO terms, added vep to bioalcidae ; https://github.com/lindenb/jvarkit/commit/96d5c5dcc556f399bb8cf34bbc4d6f31fbebc8c5
Mon Mar 13 17:51:01 2017 +0100 ; cont, textflow, thread/runner for bioalcidae ; https://github.com/lindenb/jvarkit/commit/0b8dff3e5660b68610f730a41eee38f57cd801df
Fri Apr 8 17:17:56 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/a943c4dfb102e4e4475d733fb32eb3dd22eb2760
Thu Apr 7 17:18:22 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/136836f4840412ca7268127725d3b310fc8b9f1d
Mon Feb 22 09:14:35 2016 +0100 ; bioalcidae stackoverflow ; https://github.com/lindenb/jvarkit/commit/24ed00d9f664c2c1064439b8ec888b55d939a276
Wed Feb 10 11:18:43 2016 +0100 ; flush ; https://github.com/lindenb/jvarkit/commit/4d7ef46efb7db0c610cb7a712961bf3f8160bd6e
Tue Jan 26 14:45:35 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/5201fb71870100bc4e4f40acf09d4e5f31ec4182
Fri Nov 27 15:22:25 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/04a83d5b9f0e69fd2f7087e519b0de3e2b4f9863
Thu Nov 26 12:58:31 2015 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/cfdff2e66fbeaa4627b50361a09196be2a2e1477
Fri Oct 2 16:11:52 2015 +0200 ; bioalci ; https://github.com/lindenb/jvarkit/commit/5506708d3f20dc4a6f6b9d43805f9722d47582b5
Fri Oct 2 15:13:06 2015 +0200 ; bioalci ; https://github.com/lindenb/jvarkit/commit/7408e134eeedf1a9ec7378f93ac2679ce37d57c5
Fri Jul 24 18:42:25 2015 +0200 ; bug for bam ; https://github.com/lindenb/jvarkit/commit/fe2d2d56cb59e045be890ced6c58b7afcd78d141
Mon Jul 6 16:14:07 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/ee95fe6971b5655c61d7feb22e8fa877201a9ca6
Wed May 13 15:22:57 2015 +0200 ; bioalcidae: javascript-based file (vcf...) reformatter #tweet ; https://github.com/lindenb/jvarkit/commit/0f35c1940e5f381aaf593731bb3fa8d5540fa29f
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bioalcidae** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



Bioinformatics file javascript-based reformatter ( java engine http://openjdk.java.net/projects/nashorn/ ). Something like awk for VCF, BAM, SAM, FASTQ, FASTA etc...

## About Galaxy

At first, this tool is not safe for a public Galaxy server, because the javascript code can access the filesystem.
But you can use the JVM parameter

```
-J-Djava.security.manager
```

to prevent it to access the filesystem. See [http://stackoverflow.com/questions/40177810](http://stackoverflow.com/questions/40177810)


## Variables


The program injects the following variables:


 *  out a java.io.PrintWriter ( https://docs.oracle.com/javase/7/docs/api/java/io/PrintWriter.html ) for output
 *  FILENAME a string, the name of the current input
 *   format a string, the format of the current input ("VCF"...)






#### VCF

For VCF , the program injects the following variables:

 *  header a htsjdk.variant.vcf.VCFHeader https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html
 *  iter a java.util.Iterator<htsjdk.variant.variantcontext.VariantContext>  https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html
 *  vep a com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser  https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/predictions/VepPredictionParser.java
 *  snpeff a com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser  https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/predictions/AnnPredictionParser.java





#### Fasta


 *  iter a java.util.Iterator<Fasta>





```

	public class Fasta 
		{
		public String getSequence();
		public String getName();
		public void print();
		public int getSize();
		public char charAt(int i);
		}

```






#### BLAST


 *  iter a java.util.Iterator<gov.nih.nlm.ncbi.blast.Hit> . gov.nih.nlm.ncbi.blast.Hit is defined by the Blast Document type definition (DTD). This iterator has also a method getIteration() that returns the following interface:

```
interface BlastIteration {
		public int getNum();
		public String getQueryId();
		public String getQueryDef();
		public int getQueryLen();
		}
	}
```


	





#### INSDSEQ


 *  iter a java.util.Iterator<gov.nih.nlm.ncbi.insdseq.INSDSeq> . gov.nih.nlm.ncbi.insdseq.INSDSeq is defined by the INSDSeq Document type definition (DTD).






#### XML


 *  iter a java.util.Iterator<gov.nih.nlm.ncbi.dbsnp.Rs> . gov.nih.nlm.ncbi.dbsnp.Rs is defined by the XSD schema http://ftp.ncbi.nlm.nih.gov/snp/specs/docsum_current.xsd






#### BAM or SAM



 *  header a htsjdk.samtools.SAMFileHeader http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html
 *  iter a htsjdk.samtools.SAMRecordIterator  https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecordIterator.html






#### FASTQ


 *  iter a java.util.Iterator<htsjdk.samtools.fastq.FastqRecord>  https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/fastq/FastqRecord.html






### Example



#### BAM

getting an histogram of the length of the reads


```

L={};
while(iter.hasNext()) {
	var rec=iter.next();
	if( rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) continue;
	var n= rec.getReadLength();
	if(n in L) {L[n]++;} else {L[n]=1;}
	}
for(var i in L) {
	out.println(""+i+"\t"+L[i]);
	}

```





#### Fasta

"Creating a consensus based on 'x' number of fasta files" ( https://www.biostars.org/p/185162/#185168)



```

$ echo -e ">A_2795\nTCAGAAAGAACCTC\n>B_10\nTCAGAAAGCACCTC\n>C_3\nTCTGAAAGCACTTC" |\
java -jar ~/src/jvarkit-git/dist/bioalcidae.jar -F fasta -e 'var consensus=[];while(iter.hasNext()) { var seq=iter.next();out.printlnseq.name+"\t"+seq.sequence);for(var i=0;i< seq.length();++i) {while(consensus.length <= i) consensus.push({}); var b = seq.charAt(i);if(b in consensus[i]) {consensus[i][b]++;} else {consensus[i][b]=1;} } } out.print("Cons.\t"); for(var i=0;i< consensus.length;i++) {var best=0,base="N"; for(var b in consensus[i]) {if(consensus[i][b]>best) { best= consensus[i][b];base=b;}} out.print(base);} out.println();' 


A_2795      TCAGAAAGAACCTC
B_10         TCAGAAAGCACCTC
C_3           TCTGAAAGCACTTC
Cons.        TCAGAAAGCACCTC

```





#### VCF

Reformating a VCF
we want to reformat a VCF with header


```

CHROM POS REF ALT GENOTYPE_SAMPLE1 GENOTYPE_SAMPLE2 ... GENOTYPE_SAMPLEN

```


we use the following javascript file:



```

var samples = header.sampleNamesInOrder;
out.print("CHROM\tPOS\tREF\tALT");
for(var i=0;i< samples.size();++i)
	{
	out.print("\t"+samples.get(i));
	}
out.println();

while(iter.hasNext())
	{
	var ctx = iter.next();
	if(ctx.alternateAlleles.size()!=1) continue;
	out.print(ctx.getContig() +"\t"+ctx.start+"\t"+ctx.reference.displayString+"\t"+ctx.alternateAlleles.get(0).displayString);
	for(var i=0;i< samples.size();++i)
		{
		var g = ctx.getGenotype(samples.get(i));

		out.print("\t");

		if(g.isHomRef())
			{
			out.print("0");
			}
		else if(g.isHomVar())
			{
			out.print("2");
			}
		else if(g.isHet())
			{
			out.print("1");
			}
		else
			{
			out.print("-9");
			}
		}
	out.println();
	}

```






```

$ curl -s  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" | \
gunzip -c | java -jar ./dist/bioalcidae.jar -f jeter.js -F vcf | head -n 5 | cut -f 1-10

CHROM	POS	REF	ALT	HG00096	HG00097	HG00099	HG00100	HG00101	HG00102
22	16050075	A	G	0	0	0	0	0	0
22	16050115	G	A	0	0	0	0	0	0
22	16050213	C	T	0	0	0	0	0	0
22	16050319	C	T	0	0	0	0	0	0

```



***
for 1000 genome data, print CHROM/POS/REF/ALT/AF(europe):



```

$ curl  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz" |  gunzip -c |\
java -jar dist/bioalcidae.jar  -F VCF -e 'while(iter.hasNext()) {var ctx=iter.next(); if(!ctx.hasAttribute("EUR_AF") || ctx.alternateAlleles.size()!=1) continue; out.println(ctx.getContig()+"\t"+ctx.start+"\t"+ctx.reference.displayString+"\t"+ctx.alternateAlleles.get(0).displayString+"\t"+ctx.getAttribute("EUR_AF"));}' 

1	10177	A	AC	0.4056
1	10235	T	TA	0
1	10352	T	TA	0.4264
1	10505	A	T	0
1	10506	C	G	0
1	10511	G	A	0
1	10539	C	A	0.001
1	10542	C	T	0
1	10579	C	A	0
1	10616	CCGCCGTTGCAAAGGCGCGCCG	C	0.994
(...)

```





#### Blast

```
$ cat ~/input.blastn.xml | java -jar dist/bioalcidae.jar -F blast -e 'while(iter.hasNext())
 	{
 	var query  = iter.getIteration();
 	var hit = iter.next();
 	out.println(query.getQueryDef()+" Hit: "+hit.getHitDef()+"  num-hsp = "+hit.getHitHsps().getHsp().size());
 	}'

```

output:


```
$
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Escherichia coli genome assembly FHI90 ,scaffold scaffold-6_contig-25.0_1_5253_[organism:Escherichia  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Acinetobacter baumannii AC12, complete genome  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Escherichia coli genome assembly FHI7 ,scaffold scaffold-5_contig-23.0_1_5172_[organism:Escherichia  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Escherichia coli genome assembly FHI92 ,scaffold scaffold-6_contig-18.0_1_5295_[organism:Escherichia  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Amycolatopsis lurida NRRL 2430, complete genome  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Escherichia coli genome assembly FHI87 ,scaffold scaffold-4_contig-19.0_1_5337_[organism:Escherichia  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Desulfitobacterium hafniense genome assembly assembly_v1 ,scaffold scaffold9  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Escherichia coli genome assembly FHI79 ,scaffold scaffold-4_contig-23.0_1_3071_[organism:Escherichia  num-hsp = 1
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Escherichia coli genome assembly FHI24 ,scaffold scaffold-8_contig-33.0_1_3324_[organism:Escherichia  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Escherichia coli genome assembly FHI89 ,scaffold scaffold-8_contig-14.0_1_3588_[organism:Escherichia  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Sphingobacterium sp. PM2-P1-29 genome assembly Sequencing method ,scaffold BN1088_Contig_19  num-hsp = 2
Enterobacteria phage phiX174 sensu lato, complete genome Hit: Escherichia coli genome assembly FHI43 ,scaffold scaffold-3_contig-14.0_1_2537_[organism:Escherichia  num-hsp = 1
(...)
```





#### NCBI Sequence INDSeq

```
$  curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=25,26,27&rettype=gbc" |\
java -jar dist/bioalcidae.jar  -F INSDSEQ -e 'while(iter.hasNext()) {var seq= iter.next(); out.println(seq.getINSDSeqDefinition()+" LENGTH="+seq.getINSDSeqLength());}'
```


output:

```
Blue Whale heavy satellite DNA LENGTH=422
Blue Whale heavy satellite DNA LENGTH=416
B.physalus gene for large subunit rRNA LENGTH=518
```


#### NCBI DBSNP



```
$  curl -s "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/XML/ds_chMT.xml.gz" | gunzip -c |\
 java -jar dist/bioalcidae.jar  -F dbsnp -e 'while(iter.hasNext())
 	{ var rs= iter.next();
 	out.println("rs"+rs.getRsId()+" "+rs.getSnpClass()+" "+rs.getMolType());
 	}'

rs8936 snp genomic
rs9743 snp genomic
rs1015433 snp genomic
rs1029272 snp genomic
rs1029293 snp genomic
rs1029294 snp genomic
rs1041840 snp genomic
rs1041870 snp genomic
rs1064597 snp cDNA
rs1116904 snp genomic
(...)

```




