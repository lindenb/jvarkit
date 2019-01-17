/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.bioalcidae;


import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.LineIterator;

import java.io.ByteArrayInputStream;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.Arrays;
import java.util.List;

import javax.script.Bindings;
import javax.script.CompiledScript;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.google.gson.JsonElement;
import com.google.gson.JsonNull;
import com.google.gson.JsonParser;

import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.dbsnp.Rs;
import gov.nih.nlm.ncbi.insdseq.INSDSeq;

/**

BEGIN_DOC


Bioinformatics file javascript-based reformatter ( java engine http://openjdk.java.net/projects/nashorn/ ). Something like awk for VCF, BAM, SAM, FASTQ, FASTA etc...


## Why  this name, 'BioAlcidae' ?

As 'bioalcidae' looks like an 'awk' for bioinformatics, we used '[Alcidae](https://en.wikipedia.org/wiki/Alcidae)', the taxonomic Family of the '[auk](https://en.wikipedia.org/wiki/Auk)' species.

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



END_DOC
*/


@Program(name="bioalcidae",
	description="javascript version of awk for bioinformatics",
	keywords={"sam","bam","vcf","javascript","js","nashorn"},
	biostars={276219,257346,183197,185162,153060,152016,152720,152820,218444,224402,241751,240452,248385,186610,242127,167389,187494,183197,152820,178004,156250,202400,183982,173201},
	references="\"bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files\" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734)."
	)
public class BioAlcidae
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BioAlcidae.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-F","--format"},description="force format: one of VCF BAM SAM FASTQ FASTA BLAST DBSNP. BLAST is BLAST XML version 1. DBSNP is XML output of NCBI dbsnp. INSDSEQ is XML output of NCBI EFetch rettype=gbc.")
	private String formatString = null;

	@Parameter(names={"-J","--json"},description="Optional. Reads a JSON File using google gson (https://google-gson.googlecode.com/svn/trunk/gson/docs/javadocs/index.html ) and injects it as 'userData' in the javascript context.")
	private File jsonFile = null;
	
    @Parameter(names={"-f","--scriptfile"},description="Javascript file")
    private File javascriptFile=null;
    @Parameter(names={"-e","--expression"},description="Javascript expression")
    private String javascriptExpr=null;

	
	@SuppressWarnings("unused")
	private static gov.nih.nlm.ncbi.blast.ObjectFactory _fooljavac1 = null;
	@SuppressWarnings("unused")
	private static gov.nih.nlm.ncbi.insdseq.ObjectFactory _fooljavac2 = null;
	@SuppressWarnings("unused")
	private static gov.nih.nlm.ncbi.dbsnp.ObjectFactory _fooljavac3 = null;

	private enum FORMAT {
		VCF{
			@Override
			boolean canAs(final String src) {
				return src!=null && (Arrays.asList(IOUtil.VCF_EXTENSIONS).stream().anyMatch(EXT->src.endsWith(EXT)) );
			}
			},
		SAM{
			@Override
			boolean canAs(String src) {
				return src!=null && (src.endsWith(".sam") || src.endsWith(".bam") );
			}
			},
		BAM{
			@Override
			boolean canAs(String src) {
				return src!=null && (src.endsWith(".sam") || src.endsWith(".bam") );
			}
			},
		FASTA{
			@Override
			boolean canAs(String src) {
				return src!=null && (src.endsWith(".fa") || src.endsWith(".fasta") || src.endsWith(".fa.gz") || src.endsWith(".fasta.gz")  );
			}
			},
		FASTQ{
				@Override
				boolean canAs(String src) {
					return src!=null && (src.endsWith(".fq") || src.endsWith(".fastq") || src.endsWith(".fq.gz") || src.endsWith(".fastq.gz")  );
				}
				}
			,
		BLAST{
				@Override
				boolean canAs(String src) {
					return src!=null && (src.endsWith(".blast.xml")  );
				}
				},
		INSDSEQ{
				@Override
				boolean canAs(String src) {
					return src!=null && (src.endsWith(".insdseq.xml")  );
				}
				},
		DBSNP {
					@Override
					boolean canAs(String src) {
						return src!=null && (src.endsWith(".dbsnp.xml")  );
					}
				}
			;
		abstract boolean canAs(String src);
		};
		
	
		
		private Bindings bindings=null;
		private CompiledScript  script=null;
		private FORMAT format= null;
		private PrintWriter writer=null;


		

	
	/** moved to public for knime */
	public  int executeAsVcf(final String source) throws IOException
		{
		LOG.info("source: "+source);
		VCFIterator in=null;
		try {
			in = VCFUtils.createVCFIterator(source);
			return executeAsVcf(in);
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	public  int executeAsVcf(final VCFIterator in) throws IOException
		{
		try {
			bindings.put("codec",VCFUtils.createDefaultVCFCodec());
			bindings.put("header",in.getHeader());
			bindings.put("iter",in);
			bindings.put("format","vcf");
			bindings.put("snpeff",new AnnPredictionParserFactory(in.getHeader()).get());
			bindings.put("vep",new VepPredictionParserFactory(in.getHeader()).get());
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (final Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			bindings.remove("format");
			bindings.remove("codec");
			bindings.remove("header");
			bindings.remove("iter");
			}
		}

	
	
	private int execute_bam(String source) throws IOException
		{
		SamReader in=null;
		SAMRecordIterator iter=null;
		try {
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			if(source==null)
				{
				in= srf.open(SamInputResource.of(stdin()));
				}
			else
				{
				in= srf.open(SamInputResource.of(source));
				}
			iter = in.iterator();
			bindings.put("header",in.getFileHeader());
			bindings.put("iter",iter);
			bindings.put("format","sam");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(iter);
			bindings.remove("header");
			bindings.remove("iter");
			bindings.remove("format");
			}
		}
	
	public class Fasta 
		{
		private String sequence;
		private String name;
		public String getSequence() {
			return sequence;
			}
		public String getName() {
			return name;
			}
		public void print()
			{
			BioAlcidae.this.writer.print(">");
			BioAlcidae.this.writer.print(name);
			for(int i=0;i< sequence.length();++i)
				{
				if(i%60==0) BioAlcidae.this.writer.println();
				BioAlcidae.this.writer.print(this.charAt(i));
				}
			BioAlcidae.this.writer.println();
			}
		public int length()
			{
			return sequence.length();
			}
		public int getLength()
			{
			return sequence.length();
			}
		public int size()
			{
			return sequence.length();
			}
		public int getSize()
			{
			return sequence.length();
			}
		public char charAt(int i)
			{
			return sequence.charAt(i);
			}
		@Override
		public String toString() {
			return sequence;
			}
		}
	
	public class FastaIterator extends AbstractIterator<Fasta>
		{
		private LineIterator in=null;
		@Override
		protected Fasta advance()
			{
			Fasta f=null;
			for(;;)
				{
				if(!in.hasNext()) return null;
				String line=in.next();
				if(line.trim().isEmpty()) continue;
				if(!line.startsWith(">")) throw new RuntimeException("Expected '>'. Bad fasta line :"+line);
				f=new Fasta();
				f.name=line.substring(1);
				break;
				}
			StringBuilder sb=new StringBuilder();
			while(in.hasNext())
				{
				final String line=in.peek();
				if(line.trim().isEmpty()) {in.next();continue;}
				if(line.startsWith(">")) break;//new sequence
				sb.append(in.next().trim());
				}
			f.sequence = sb.toString();
			return f;
			}
		}
	
	private  int execute_fasta(String source) throws IOException
		{
		FastaIterator iter=new FastaIterator();
		try {
			iter.in = (source==null?IOUtils.openStreamForLineIterator(stdin()):IOUtils.openURIForLineIterator(source));
			bindings.put("iter",iter);
			bindings.put("format","fasta");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter.in);
			bindings.remove("iter");
			bindings.remove("format");
			}
		}
	
	private  int execute_fastq(String source) throws IOException
		{
		InputStream in=null;
		FourLinesFastqReader iter=null;
		try {
			if(source==null)
				{
				in= stdin();
				}
			else
				{
				in= IOUtils.openURIForReading(source);
				}
			iter = new FourLinesFastqReader(in);
			bindings.put("iter",in);
			bindings.put("format","fastq");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(in);
			bindings.remove("header");
			bindings.remove("iter");
			bindings.remove("format");
			}
		}
	
	public static class BlastIteration {
		private int num=0;
		private String queryId=null;
		private String queryDef=null;
		private int queryLen=-1;
		public int getNum() {
			return num;
		}
		public String getQueryId() {
			return queryId;
		}
		public String getQueryDef() {
			return queryDef;
		}
		public int getQueryLen() {
			return queryLen;
		}
	}
	
	
	public static abstract class AbstractXMLIterator<T>  extends AbstractIterator<T> implements Closeable {
		protected final InputStream in;
		protected final JAXBContext jc;
		protected final Unmarshaller unmarshaller;
		protected XMLEventReader r=null;
		AbstractXMLIterator(final InputStream in,final String jaxbPath,final Boolean namespaceaware) throws JAXBException,XMLStreamException  {
			this.in = in;
			this.jc = JAXBContext.newInstance(jaxbPath);
			this.unmarshaller =jc.createUnmarshaller();
			final XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE,namespaceaware);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			xmlInputFactory.setXMLResolver(new XMLResolver()
				{
				@Override
				public Object resolveEntity(String publicID,
						String systemID, String baseURI, String namespace)
						throws XMLStreamException {
						LOG.warn("ignoring resolveEntity "+systemID+" "+baseURI+" "+namespace);
							return new ByteArrayInputStream(new byte[0]);
						}
				});
			final StreamSource streamSource = new StreamSource(in);
			r=xmlInputFactory.createXMLEventReader(streamSource);
			}
		
		@Override
		public void close() throws IOException {
			CloserUtil.close(r);
			CloserUtil.close(in);
			}
		protected T simpleAdvance(final Class<T> clazz,final String tagName) {
			try {
				while(r.hasNext())
					{
					final XMLEvent evt=r.peek();
					if(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals(tagName))
						{
						return unmarshaller.unmarshal(r,clazz).getValue();
						}
					else
						{
						r.next();//consumme
						}
					}
				return null;
			} catch (XMLStreamException|JAXBException e) {
				throw new RuntimeIOException(e);
			}
			
		}
	}
	
	public static class BlastIterator extends AbstractXMLIterator<Hit> {
		private BlastIteration iteration= new BlastIteration();
		BlastIterator(final InputStream in) throws JAXBException,XMLStreamException {
			super(in,"gov.nih.nlm.ncbi.blast",Boolean.TRUE);
			}
		
		public BlastIteration getIteration() {
			return this.iteration;
		}
		
		@Override
		protected Hit advance() {
			try {
				while(r.hasNext())
					{
					XMLEvent evt=r.peek();
					if(evt.isStartElement() )
						{
						final String localName= evt.asStartElement().getName().getLocalPart();
						if(localName.equals("Hit"))
							{
							return unmarshaller.unmarshal(r,Hit.class).getValue();
							}
						else if(localName.equals("Iteration"))
							{
							this.iteration= new BlastIteration();
							r.next();//consumme
							}
						else if(localName.equals("Iteration_iter-num"))
							{
							r.next();//consumme
							this.iteration.num= Integer.parseInt(r.getElementText());
							}
						else if(localName.equals("Iteration_query-ID"))
							{
							r.next();//consumme
							this.iteration.queryId= r.getElementText();
							}
						else if(localName.equals("Iteration_query-def"))
							{
							r.next();//consumme
							this.iteration.queryDef= r.getElementText();
							}
						else if(localName.equals("Iteration_query-len"))
							{
							r.next();//consumme
							this.iteration.queryLen = Integer.parseInt(r.getElementText());
							}
						else
							{
							r.next();//consumme
							}
						}
					else
						{
						r.next();//consumme
						}
					}
				return null;
			} catch (XMLStreamException|JAXBException e) {
				throw new RuntimeIOException(e);
			}
			}
		
		}
	
	private  int execute_blast(String source) throws IOException
		{
		BlastIterator iter=null;
		try {
			iter =new BlastIterator(source==null?stdin():IOUtils.openURIForReading(source));
			bindings.put("iter",iter);
			bindings.put("format","blast");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			bindings.remove("iter");
			bindings.remove("format");
			}
		}

	
	public static class InsdSeqIterator extends AbstractXMLIterator<INSDSeq>   {
		InsdSeqIterator(final InputStream in) throws JAXBException,XMLStreamException {
		super(in,"gov.nih.nlm.ncbi.insdseq",Boolean.TRUE);
		}
	
			@Override
			protected INSDSeq advance() { return super.simpleAdvance(INSDSeq.class, "INSDSeq"); 
			}
		}

	
	
	private  int execute_insdseq(String source) throws IOException
		{
		InsdSeqIterator iter=null;
		try {
			iter =new InsdSeqIterator(source==null?stdin():IOUtils.openURIForReading(source));
			bindings.put("iter",iter);
			bindings.put("format","insdseq");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			bindings.remove("iter");
			bindings.remove("format");
			}
		}

	
	private  int execute_dbsnp(String source) throws IOException
		{
		AbstractXMLIterator<Rs> iter=null;
		try {
			iter =new AbstractXMLIterator<Rs>(
					source==null?stdin():IOUtils.openURIForReading(source),
					"gov.nih.nlm.ncbi.dbsnp",
					Boolean.TRUE)
				{
				@Override
				protected Rs advance() {
					return super.simpleAdvance(Rs.class, "Rs");
					}
				};
			bindings.put("iter",iter);
			bindings.put("format","dbsnp");
			this.script.eval(bindings);
			this.writer.flush();
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			bindings.remove("iter");
			bindings.remove("format");
			}
		}

	
	
	private  int execute(FORMAT fmt,String source) throws IOException
		{
		switch(fmt)
			{
			case VCF: return executeAsVcf(source);
			case BAM: case SAM: return execute_bam(source);
			case FASTQ: return execute_fastq(source);
			case FASTA: return execute_fasta(source);
			case BLAST: return execute_blast(source);
			case INSDSEQ: return execute_insdseq(source);
			case DBSNP: return execute_dbsnp(source);
			default: throw new IllegalStateException("Bad file format "+fmt);
			}
		}

	
	@Override
	public int doWork(final List<String> args) {
		
		Reader jsonIn = null;
		
		if(this.formatString!=null)
			{
			try {
				this.format=FORMAT.valueOf(this.formatString.toUpperCase());
				} catch (Exception err) {
				LOG.error(err);
				return -1;
				}
			}
		
		try
			{
			this.script = super.compileJavascript(this.javascriptExpr,this.javascriptFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		
		try
			{
			
			this.bindings = this.script.getEngine().createBindings();
			this.writer = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			this.bindings.put("out",this.writer);
			
			if( this.jsonFile != null) {
				LOG.info("Reading JSON FILE "+this.jsonFile);
				jsonIn = IOUtils.openFileForBufferedReading(this.jsonFile);
				final JsonParser gsonParser = new JsonParser();
				final JsonElement jsonElement = gsonParser.parse(jsonIn);
				jsonIn.close();
				jsonIn=null;
				this.bindings.put("userData", jsonElement);
			} else
				{
				this.bindings.put("userData", JsonNull.INSTANCE);
				}
			
			if(args.isEmpty()  )
				{
				if(this.format==null)
					{
					LOG.error("Format must be specified when reading from stdin");
					return -1;
					}
				return execute(this.format,null);
				}
			for(final String filename:IOUtils.unrollFiles(args))
				{
				this.bindings.put("FILENAME",filename);
				FORMAT fmt=this.format;
				if(fmt==null)
					{
					for(FORMAT t:FORMAT.values())
						{
						if(t.canAs(filename))
							{
							fmt=t;
							break;
							}
						}
					}
				if(fmt==null)
					{
					LOG.error("Cannot get file format for "+filename+". Please specifiy using option -F");
					return -1;
					}
				final int errors= execute(fmt,filename);
				if(errors!=0)
					{
					return errors;
					}
				}
				
			
			this.writer.flush();
			if(this.writer.checkError())
				{
				LOG.warn("I/O errror. Something wrong happended");
				}
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(this.writer!=null) this.writer.flush();
			if(this.outputFile!=null)
				{
				CloserUtil.close(this.writer);
				}
			this.writer=null;
			this.bindings=null;
			this.bindings=null;
			this.format=null;
			CloserUtil.close(jsonIn);
			}
		}
	
	
	
	
	
	
	
	public static void main(String[] args) {
		new BioAlcidae().instanceMainWithExit(args);
	}
}
