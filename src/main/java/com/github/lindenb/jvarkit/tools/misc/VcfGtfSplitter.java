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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

/**

BEGIN_DOC

### Example

```
$ java -jar dist/vcfgtfsplitter.jar  -m jeter.manifest --gtf  input.gtf.gz -o jeter.zip src/test/resources/test_vcf01.vcf 

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     4741  2019-11-18 10:52   d2/a0884d9bf86378a2cd0bfacef19723/ENSG00000188157.vcf.gz
     4332  2019-11-18 10:52   4d/11ff4f2413d4a369ee9a51192acad3/ENSG00000131591.vcf.gz
---------                     -------
     9073                     2 files

$ column -t jeter.manifest 
#chrom  start    end      Gene-Id          Gene-Name  Gene-Biotype    path                                                      Count_Variants
1       955502   991496   ENSG00000188157  AGRN       protein_coding  d2/a0884d9bf86378a2cd0bfacef19723/ENSG00000188157.vcf.gz  20
1       1017197  1051741  ENSG00000131591  C1orf159   protein_coding  4d/11ff4f2413d4a369ee9a51192acad3/ENSG00000131591.vcf.gz  6


$ java -jar dist/vcfgtfsplitter.jar -T --index  --gtf  jeter.gtf  -m jeter.manifest -o jeter.zip src/test/resources/test_vcf01.vcf
$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     4749  2019-11-18 10:56   93/00d3aa560bdbe3b3ff964e24303107/ENST00000379370.vcf.gz
     4749  2019-11-18 10:56   93/00d3aa560bdbe3b3ff964e24303107/ENST00000379370.vcf.gz.tbi
     3108  2019-11-18 10:56   6d/001c37d11b192bcd77ac089ec258f2/ENST00000477585.vcf.gz
(...)
     2969  2019-11-18 10:56   3a/8e697b04e03716942b362afbe7ee51/ENST00000472741.vcf.gz
     2969  2019-11-18 10:56   3a/8e697b04e03716942b362afbe7ee51/ENST00000472741.vcf.gz.tbi
     2969  2019-11-18 10:56   2c/4ced556d7722b978453c38d5525ee0/ENST00000480643.vcf.gz
     2969  2019-11-18 10:56   2c/4ced556d7722b978453c38d5525ee0/ENST00000480643.vcf.gz.tbi
---------                     -------
   175030                     46 files

$ column -t jeter.manifest 
#chrom  start    end      Gene-Id          Gene-Name  Gene-Biotype    Transcript-Id    path                                                      Count_Variants
1       955502   991496   ENSG00000188157  AGRN       protein_coding  ENST00000379370  93/00d3aa560bdbe3b3ff964e24303107/ENST00000379370.vcf.gz  20
1       969485   976105   ENSG00000188157  AGRN       protein_coding  ENST00000477585  6d/001c37d11b192bcd77ac089ec258f2/ENST00000477585.vcf.gz  3
1       970246   976777   ENSG00000188157  AGRN       protein_coding  ENST00000469403  a9/fe6f84e1a425797d5b58ceae3a8313/ENST00000469403.vcf.gz  2
1       983908   984774   ENSG00000188157  AGRN       protein_coding  ENST00000492947  55/b5a514c94417efe2632aba379d8247/ENST00000492947.vcf.gz  1
1       1017197  1051461  ENSG00000131591  C1orf159   protein_coding  ENST00000379339  7b/26f9faff411f6e056fefc89cd15a24/ENST00000379339.vcf.gz  6
1       1017197  1051736  ENSG00000131591  C1orf159   protein_coding  ENST00000448924  56/fad8194b90771b0c3e931bc4772be1/ENST00000448924.vcf.gz  6
1       1017197  1051736  ENSG00000131591  C1orf159   protein_coding  ENST00000294576  df/3aa6bebac650a6c9f59409521ae17d/ENST00000294576.vcf.gz  6
(...)
```

END_DOC
*/
@Program(
		name="vcfgtfsplitter",
		description="Split VCF+VEP by gene/transcript using a GTF file.",
		creationDate="20191118",
		modificationDate="20191118",
		keywords= {"genes","vcf","split","gtf"}
		)
public class VcfGtfSplitter
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfGtfSplitter.class).make();
	
	
	@Parameter(names={"-o","--output"},description= ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-m","--manifest"},description="Manifest Bed file output containing chrom/start/end of each gene")
	private Path manifestFile = null;
	@Parameter(names={"-T","--transcript"},description= "split by transcript. (default is to split per gene)")
	private boolean split_by_transcript = false;
	@Parameter(names={"-g","-G","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfPath = null;
	@Parameter(names={"--ignore-filtered"},description="Ignore FILTERED variant")
	private boolean ignoreFiltered = false;
	@Parameter(names={"-C","--contig","--chromosome"},description="Limit to those contigs.")
	private Set<String> limitToContigs = new HashSet<>();
	@Parameter(names={"--bcf"},description="Use bcf format")
	private boolean use_bcf= false;
	@Parameter(names={"--index"},description="index files")
	private boolean index_vcf= false;
	@Parameter(names={"--features"},description="Features to keep. Comma separated. a set of 'cds,exon,intron,transcript,utr,utr5,utr3,stop,start'")
	private String featuresString = "cds,exon,intron,transcript,utr5,utr3,stop,start";


	
	private boolean use_cds = false;
	private boolean use_intron = false;
	private boolean use_exon = false;
	private boolean use_stop = false;
	private boolean use_start = false;
	private boolean use_utr5 = false;
	private boolean use_utr3 = false;

	
	private abstract class AbstractSplitter
		{
		abstract void addMetadata(VCFHeader h);
		abstract Locatable getInterval();
		abstract boolean accept(final VariantContext ctx);
		abstract String getId();
		abstract void printManifest(PrintWriter pw);
		}
	private  class TranscriptSplitter extends AbstractSplitter
		{
		final Transcript transcript;
		TranscriptSplitter(final Transcript transcript) {
			this.transcript = transcript;
			}
		@Override
		Locatable getInterval() {
			return this.transcript;
			}
		@Override
		void addMetadata(VCFHeader h) {
			final Gene gene = transcript.getGene();
			h.addMetaDataLine(new VCFHeaderLine("split.gene-id", gene.getId()));
			h.addMetaDataLine(new VCFHeaderLine("split.gene-name", gene.getGeneName()));
			h.addMetaDataLine(new VCFHeaderLine("split.gene-biotype", gene.getGeneBiotype()));
			h.addMetaDataLine(new VCFHeaderLine("split.transcript-id", transcript.getId()));
			}
		@Override
		boolean accept(VariantContext ctx) {
			return testTranscript(this.transcript,ctx);
			}
		@Override
		String getId() {
			return transcript.getId();
			}
		@Override
		void printManifest(PrintWriter pw) {
			final Gene gene = transcript.getGene();
			pw.print("\t");
			pw.print(gene.getId());
			pw.print("\t");
			pw.print(gene.getGeneName());
			pw.print("\t");
			pw.print(gene.getGeneBiotype());
			pw.print("\t");
			pw.print(transcript.getId());
			}
		}
	private  class GeneSplitter extends AbstractSplitter
		{
		final Gene gene;
		GeneSplitter(final Gene gene) {
			this.gene = gene;
			}
		@Override
		Locatable getInterval() {
			return this.gene;
			}
		@Override
		void addMetadata(VCFHeader h) {
			h.addMetaDataLine(new VCFHeaderLine("split.gene-id", gene.getId()));
			h.addMetaDataLine(new VCFHeaderLine("split.gene-name", gene.getGeneName()));
			h.addMetaDataLine(new VCFHeaderLine("split.gene-biotype", gene.getGeneBiotype()));
			}
		@Override
		boolean accept(final VariantContext ctx) {
			return gene.getTranscripts().stream().anyMatch(T->testTranscript(T, ctx));
			}
		@Override
		String getId() {
			return gene.getId();
			}
		@Override
		void printManifest(PrintWriter pw) {
			pw.print("\t");
			pw.print(gene.getId());
			pw.print("\t");
			pw.print(gene.getGeneName());
			pw.print("\t");
			pw.print(gene.getGeneBiotype());
			}
		}
		
	public VcfGtfSplitter()
		{
		
		}
	
	private boolean testTranscript(final Transcript transcript,final VariantContext ctx) {
		if(this.use_exon && transcript.hasExon() && transcript.getExons().stream().anyMatch(FEAT->FEAT.overlaps(ctx))) return true;
		if(this.use_intron && transcript.hasIntron() && transcript.getIntrons().stream().anyMatch(FEAT->FEAT.overlaps(ctx))) return true;
		if(this.use_cds && transcript.hasCDS() && transcript.getAllCds().stream().anyMatch(FEAT->FEAT.overlaps(ctx))) return true;
		if(this.use_utr5) {
			if(transcript.isPositiveStrand() && transcript.getUTR5().isPresent() && transcript.getUTR5().get().overlaps(ctx)) return true;
			if(transcript.isNegativeStrand() && transcript.getUTR3().isPresent() && transcript.getUTR3().get().overlaps(ctx)) return true;
			}
		if(this.use_utr3) {
			if(transcript.isPositiveStrand() && transcript.getUTR3().isPresent() && transcript.getUTR3().get().overlaps(ctx)) return true;
			if(transcript.isNegativeStrand() && transcript.getUTR5().isPresent() && transcript.getUTR5().get().overlaps(ctx)) return true;
			}
		if(this.use_stop && transcript.hasCodonStopDefined() && transcript.getCodonStop().get().getBlocks().stream().anyMatch(FEAT->FEAT.overlaps(ctx))) return true;
		if(this.use_start && transcript.hasCodonStartDefined() && transcript.getCodonStart().get().getBlocks().stream().anyMatch(FEAT->FEAT.overlaps(ctx))) return true;
		
		return false;
	}
	
	private void split(
			final AbstractSplitter splitter,
			final VCFFileReader vcfFileReader,
			final VCFHeader header,
			final SAMSequenceDictionary dict,
			final ArchiveFactory archiveFactory,
			final Path tmpVcf,
			final PrintWriter manifest
			) throws IOException {
		
		final Locatable interval = splitter.getInterval();
		CloseableIterator<VariantContext> iter=vcfFileReader.query(interval);
		if(!iter.hasNext()) {
			iter.close();
			return;
		}
		
		
		final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
		vcwb.setCreateMD5(false);
		vcwb.setReferenceDictionary(null);
		final TabixIndexCreator tabixIndexCreator;
		final Path tbiPath;
		if(this.index_vcf) {
			if(dict==null) {
				tabixIndexCreator = new TabixIndexCreator(TabixFormat.VCF);
			} else {
				tabixIndexCreator = new TabixIndexCreator(dict, TabixFormat.VCF);
			}
			vcwb.setIndexCreator(tabixIndexCreator);
			vcwb.setOption(Options.INDEX_ON_THE_FLY);
			tbiPath =tmpVcf.getParent().resolve(tmpVcf.getFileName().toString()+FileExtensions.TABIX_INDEX);
			} 
		else {
			tabixIndexCreator = null;
			tbiPath = null;
			}
		
		vcwb.setOutputPath(tmpVcf);
		
		
		final VariantContextWriter out = vcwb.build();
		final VCFHeader header2=new VCFHeader(header);
		super.addMetaData(header2);
		splitter.addMetadata(header2);
		out.writeHeader(header2);

		
		int count_ctx = 0;
		while(iter.hasNext()) {
			final VariantContext ctx = iter.next();
			if(this.ignoreFiltered && ctx.isFiltered()) continue;
			if(!splitter.accept(ctx)) continue;
			count_ctx++;
			out.add(ctx);
			}
		iter.close();
		out.close();

		if(count_ctx==0) {
			Files.delete(tmpVcf);
			if(tbiPath!=null) Files.delete(tbiPath);
			return ;
			}
		
		final String md5 = StringUtils.md5(interval.getContig()+":"+splitter.getId());
		final String filename =  md5.substring(0,2) + File.separatorChar + md5.substring(2) + File.separator+splitter.getId().replaceAll("[/\\:]", "_") + (this.use_bcf?FileExtensions.BCF:FileExtensions.COMPRESSED_VCF);

		
		OutputStream os = archiveFactory.openOuputStream(filename);
		IOUtils.copyTo(tmpVcf, os);
		os.flush();
		os.close();
		if(tbiPath!=null) {
			os = archiveFactory.openOuputStream(filename+FileExtensions.TABIX_INDEX);
			IOUtils.copyTo(tmpVcf, os);
			os.flush();
			os.close();
			Files.delete(tbiPath);
			}
		Files.delete(tmpVcf);
		
		
		manifest.print(interval.getContig());
		manifest.print('\t');
		manifest.print(interval.getStart()-1);
		manifest.print('\t');
		manifest.print(interval.getEnd());
		splitter.printManifest(manifest);
		manifest.print('\t');
		manifest.print((archiveFactory.isZip()?"":this.outputFile.toString()+File.separator)+filename);
		manifest.print('\t');
		manifest.println(count_ctx);
		}
	@Override
	public int doWork(final List<String> args) {
		ArchiveFactory archiveFactory = null;
		PrintWriter manifest = null;
		VCFFileReader vcfFileReader = null;
		try {
			for(final String s: featuresString.split("[;, ]")) {
				if(StringUtils.isBlank(s)) continue;
				if(s.equals("cds")) { use_cds = true; }
				else if(s.equals("intron")) { use_cds = true; }
				else if(s.equals("exon")) { use_exon = true; }
				else if(s.equals("stop")) { use_stop = true; }
				else if(s.equals("start")) { use_start = true; }
				else if(s.equals("transcript")) { use_exon = true; use_intron = true; }
				else if(s.equals("utr5")) { use_utr5= true;}
				else if(s.equals("utr3")) { use_utr3= true;}
				else if(s.equals("utr")) { use_utr3= true;use_utr5= true;}
				else {
					LOG.error("unknown code "+s+" in "+this.featuresString);
					return -1;
				}
			}
			
			
			final Path tmpVcf = Files.createTempFile("tmp.",(use_bcf?FileExtensions.BCF:FileExtensions.COMPRESSED_VCF));
			String input = oneAndOnlyOneFile(args);
			vcfFileReader = new VCFFileReader(Paths.get(input),true);
			final VCFHeader header1 = vcfFileReader.getFileHeader();
			final SAMSequenceDictionary dict = header1.getSequenceDictionary();
			if(dict==null && this.use_bcf) {
				throw new JvarkitException.VcfDictionaryMissing(input);
			}
			
			if(dict!=null && !limitToContigs.isEmpty())
				{
				final ContigNameConverter ctgNameConverter = ContigNameConverter.fromOneDictionary(dict);
				final Set<String> set2 = new HashSet<>(this.limitToContigs.size());
				for(final String ctg:this.limitToContigs)  {
					final String ctg2 = ctgNameConverter.apply(ctg);
					if(StringUtils.isBlank(ctg2)) {
						LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(ctg, dict));
						return -1;
						}
					set2.add(ctg2);
					}
				this.limitToContigs = set2;
				}
			
			final List<Gene> all_genes;
			try(GtfReader gtfReader=new GtfReader(this.gtfPath)) {
				if(dict!=null) gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
				all_genes = gtfReader.
						getAllGenes().
						stream().
						filter(G->{
							if(this.limitToContigs.isEmpty()) return true;
							return this.limitToContigs.contains(G.getContig());
						}).collect(Collectors.toList());
				}
			
			archiveFactory = ArchiveFactory.open(this.outputFile);
			archiveFactory.setCompressionLevel(0);
			
			manifest = new PrintWriter(this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(manifestFile));
			manifest.print("#chrom\tstart\tend\tGene-Id\tGene-Name\tGene-Biotype");
			if(this.split_by_transcript)
				{
				manifest.print("\tTranscript-Id");
				}
			manifest.println("\tpath\tCount_Variants");

			if(this.split_by_transcript) {
				final Iterator<Transcript> triter = all_genes.
						stream().
						flatMap(G->G.getTranscripts().stream()).iterator();
				while(triter.hasNext()) {
					final Transcript tr=triter.next();
					final TranscriptSplitter splitter = new TranscriptSplitter(tr);
					split(splitter,vcfFileReader,header1,dict,archiveFactory,tmpVcf,manifest);
				}
				
			} else {
				for(Gene gene: all_genes) {
					final GeneSplitter splitter = new GeneSplitter(gene);
					split(splitter,vcfFileReader,header1,dict,archiveFactory,tmpVcf,manifest);
				}
			}
			
			vcfFileReader.close();
			vcfFileReader = null;
			manifest.flush();
			manifest.close();
			manifest = null;
			archiveFactory.close();
			Files.deleteIfExists(tmpVcf);
			return RETURN_OK;
			}
		catch(final Exception err) 
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfFileReader);
			CloserUtil.close(archiveFactory);
			CloserUtil.close(manifest);
			}
		}
	
	
		
		
	 	
	
	public static void main(final String[] args)
		{
		new VcfGtfSplitter().instanceMainWithExit(args);
		}
	}
