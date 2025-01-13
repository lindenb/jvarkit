/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfsplitgene;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC

### Example

```
java -jar dist/vcfgenesplitter.jar -o jeter.zip src/test/resources/rotavirus_rf.ann.vcf.gz -m jeter.mf

$ unzip -l jeter.zip
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     1565  2019-05-27 11:26   2c/8fb9d2539e3f30d1d9b06f9ec54c4c/Gene_18_3284.vcf.gz
     2278  2019-05-27 11:26   4e/4897c51fe2dd067a8b75c19f111477/Gene_1621_1636.vcf.gz
     2278  2019-05-27 11:26   74/ca4273c3d5803c5865891c808234da/UniProtKB_Swiss-Prot:P12472.vcf.gz
     2264  2019-05-27 11:26   23/6b59cfe4fdd33a5f4feeb55521dd34/Gene_50_2557.vcf.gz
     2169  2019-05-27 11:26   b3/4bda8d8502e64e442fce077e45ded6/Gene_9_2339.vcf.gz
     2106  2019-05-27 11:26   b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz
     2023  2019-05-27 11:26   6f/8472e9f192c92bf46e4893b2367b7e/Gene_23_1216.vcf.gz
     1862  2019-05-27 11:26   3c/513d82eaea18447dd5f621f92b40e6/Gene_0_1073.vcf.gz
     1655  2019-05-27 11:26   84/977eac8cdef861cbd3109209675d21/Gene_0_1058.vcf.gz
     1754  2019-05-27 11:26   db/aee9cc8f5c9c3d39c7af4cec63b7a5/Gene_0_1061.vcf.gz
     1746  2019-05-27 11:26   b0/133c483f0ea676f8d29ab1f2daee5d/Gene_41_568.vcf.gz
     1664  2019-05-27 11:26   59/0fd5c1e8d6d60a986a0021fe357514/Gene_20_616.vcf.gz
     1663  2019-05-27 11:26   83/bc905cf311428ab80ce59aaf503838/Gene_78_374.vcf.gz
---------                     -------
    25027                     13 files

$ cat jeter.mf
#chrom	start	end	key	path	Count_Variants
RF01	969	970	ANN/GeneId	Gene_18_3284	2c/8fb9d2539e3f30d1d9b06f9ec54c4c/Gene_18_3284.vcf.gz	1
RF02	250	1965	ANN/GeneId	Gene_1621_1636	4e/4897c51fe2dd067a8b75c19f111477/Gene_1621_1636.vcf.gz	5
RF02	250	1965	ANN/GeneId	UniProtKB/Swiss-Prot:P12472	74/ca4273c3d5803c5865891c808234da/UniProtKB_Swiss-Prot:P12472.vcf.gz	5
RF03	1220	2573	ANN/GeneId	Gene_50_2557	23/6b59cfe4fdd33a5f4feeb55521dd34/Gene_50_2557.vcf.gz	8
RF04	886	1920	ANN/GeneId	Gene_9_2339	b3/4bda8d8502e64e442fce077e45ded6/Gene_9_2339.vcf.gz	7
RF05	40	1339	ANN/GeneId	Gene_32_1507	b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz	6
RF06	516	1132	ANN/GeneId	Gene_23_1216	6f/8472e9f192c92bf46e4893b2367b7e/Gene_23_1216.vcf.gz	5
RF07	97	952	ANN/GeneId	Gene_0_1073	3c/513d82eaea18447dd5f621f92b40e6/Gene_0_1073.vcf.gz	4
RF08	925	992	ANN/GeneId	Gene_0_1058	84/977eac8cdef861cbd3109209675d21/Gene_0_1058.vcf.gz	2
RF09	293	414	ANN/GeneId	Gene_0_1061	db/aee9cc8f5c9c3d39c7af4cec63b7a5/Gene_0_1061.vcf.gz	3
RF10	45	175	ANN/GeneId	Gene_41_568	b0/133c483f0ea676f8d29ab1f2daee5d/Gene_41_568.vcf.gz	3
RF11	73	79	ANN/GeneId	Gene_20_616	59/0fd5c1e8d6d60a986a0021fe357514/Gene_20_616.vcf.gz	1
RF11	73	79	ANN/GeneId	Gene_78_374	83/bc905cf311428ab80ce59aaf503838/Gene_78_374.vcf.gz	1


```

END_DOC
*/
@Program(
		name="vcfgenesplitter",
		description="Split VCF+VEP by gene/transcript.",
		creationDate = "20160310",
		modificationDate="202220531",
		keywords= {"genes","vcf"},
		jvarkit_amalgamion =  true,
		menu="VCF Manipulation"
		)
public class VcfGeneSplitter
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfGeneSplitter.class).make();
	
	private class KeyGene{
		final String extractor;
		final String key;
		final String geneName;
		int count_variants = 0;
		Path tmpVcfPath = null;
		PrintWriter pw = null;
		VCFEncoder vcfEncoder = null;
		long lastModificationDate = 0L;
		int minPos=Integer.MAX_VALUE;
		int maxPos=0;

		
		KeyGene(final String extractor,final String key,final String gene)  throws IOException {
			this.extractor = extractor;
			this.key = key;
			this.geneName = StringUtils.isBlank(gene)?".":gene;
			}
		@Override
		public int hashCode() {
			int h = extractor.hashCode();
			h= h*31 + key.hashCode();
			//h= h*31 + gene.hashCode();
			return h;
			}
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof KeyGene)) return false;
			final KeyGene kg = KeyGene.class.cast(obj);
			return this.extractor.equals(kg.extractor) && this.key.equals(kg.key) ;
			}
	
		public void write(final VCFHeader header,final VariantContext ctx) throws IOException {
			if(this.pw!=null) {
				//nothing
				}
			else if(this.tmpVcfPath ==null ) {
				LOG.info("Opening VCF for "+this.extractor+" "+this.key+" "+this.geneName);
				this.tmpVcfPath =  Files.createTempFile("tmp.", ".vcf");
				this.pw = new PrintWriter(Files.newBufferedWriter(this.tmpVcfPath, StandardOpenOption.APPEND));
				final VCFHeader h2 = new VCFHeader(header);
				h2.addMetaDataLine(new VCFHeaderLine("GtfFileSplitter.Name",String.valueOf(this.key)));
				h2.addMetaDataLine(new VCFHeaderLine("GtfFileSplitter.Gene",String.valueOf(this.geneName)));
				JVarkitVersion.getInstance().addMetaData(VcfGeneSplitter.this, h2);

				this.vcfEncoder  = new VCFEncoder(h2, false, false);

				for(final String s: VCFUtils.convertVCFHeaderToList(h2)) {
					this.pw.println(s);
					}
				}
			else
				{
				this.pw = new PrintWriter(Files.newBufferedWriter(this.tmpVcfPath, StandardOpenOption.APPEND));
				}
			this.lastModificationDate = System.currentTimeMillis();
			this.count_variants++;
			this.minPos = Math.min(ctx.getStart(), minPos);
			this.maxPos = Math.max(ctx.getEnd(), maxPos);
			this.vcfEncoder.write(this.pw,ctx);
			this.pw.println();
			}
		@Override
		public String toString() {
			return this.extractor+" "+this.key+" "+this.geneName;
			}
		}
	
	
	
	
	
	@Parameter(names={"-o","--output"},description= ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-m","--manifest"},description="Manifest Bed file output containing chrom/start/end of each gene")
	private Path manifestFile = null;
	@Parameter(names={"-l","--list"},description= "list all available extractors", help=true)
	private boolean list_extractors = false;
	@Parameter(names={"-e","-E","--extractors"},description=GeneExtractorFactory.OPT_DESC)
	private String extractorsNames="ANN/GeneId VEP/GeneId";
	@Parameter(names={"--ignore-filtered"},description="Ignore FILTERED variant")
	private boolean ignoreFiltered = false;
	@Parameter(names={"-n","--min-variant"},description="Minimum number of variants required to write a vcf. don't write if num(variant) < 'x' ")
	private int min_number_of_ctx = 1;
	@Parameter(names={"-M","--max-variant"},description="Maximum number of variants required to write a vcf. don't write if num(variant) > 'x' . '<=0' is ignore")
	private int max_number_of_ctx = -1;
	@Parameter(names={"--open-max"},description="Maximum number of opened VCF writers at the same time.")
	private int max_open_files = 100;
		
	public VcfGeneSplitter()
		{
		
		}
	
	
	private int run(final List<String> args) {
		final List<KeyGene> keyGenes = new ArrayList<>();

		try {
			
					
		try(VCFIterator iterator = super.openVCFIterator(oneFileOrNull(args))) {
			final VCFHeader header = iterator.getHeader();
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final OrderChecker<VariantContext> order = new OrderChecker<VariantContext>(dict,false);
			final GeneExtractorFactory geneExtractorFactory = new GeneExtractorFactory(header);
			final List<GeneExtractorFactory.GeneExtractor> extractors = geneExtractorFactory.parse(this.extractorsNames);
			if(extractors.isEmpty()) {
				LOG.warn("No extractor defined!");
				return -1;
				}
			try(ArchiveFactory archiveFactory = ArchiveFactory.open(this.outputFile)) {
						try(PrintWriter manifest = new PrintWriter(this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(manifestFile))) {
							manifest.println("#chrom\tstart\tend\tsplitter\tgene\tkey\tpath\tCount_Variants");
							manifest.flush();
							String prevCtg = null;
							for(;;)
								{
								final VariantContext ctx = iterator.hasNext()?order.apply(iterator.next()):null;
								
								if(this.ignoreFiltered && ctx!=null && ctx.isFiltered()) continue;
				
								if(ctx==null || !ctx.getContig().equals(prevCtg))
									{
									for(KeyGene kg: keyGenes) {
										if(kg.pw!=null) {
											kg.pw.flush();
											kg.pw.close();
											kg.pw=null;
											}
										if ( kg.count_variants < this.min_number_of_ctx )  {
											LOG.info("skipping "+kg+" because there are not enough variants. N="+kg.count_variants+"<"+this.min_number_of_ctx);
											continue;
											}
										if ( this.max_number_of_ctx!=-1 && kg.count_variants > this.max_number_of_ctx ) {
											LOG.info("skipping "+kg+" because there are too many variants. N="+kg.count_variants+">"+this.max_number_of_ctx);
											continue;
											}
				
										final String md5 = StringUtils.md5(prevCtg+":"+kg.extractor+":"+kg.key);
										final String parentDir = md5.substring(0,2) + File.separatorChar + md5.substring(2);
										final String filename0 =
												parentDir + File.separator+
												kg.key.replaceAll("[/\\:]", "_") + ".vcf.gz";
										
										
										try(final BlockCompressedOutputStream os = new BlockCompressedOutputStream(archiveFactory.openOuputStream(filename0),(Path)null)) {
											IOUtils.copyTo(kg.tmpVcfPath, os);
											os.flush();
											}
										
										manifest.print(prevCtg);
										manifest.print('\t');
										manifest.print(kg.minPos-1);
										manifest.print('\t');
										manifest.print(kg.maxPos);
										manifest.print('\t');
										manifest.print(kg.extractor);
										manifest.print('\t');
										manifest.print(kg.geneName);
										manifest.print('\t');
										manifest.print(kg.key);
										manifest.print('\t');
										manifest.print(
											archiveFactory.isTarOrZipArchive()?
											filename0:
											this.outputFile.resolve(filename0).toAbsolutePath().toString()
											);
										manifest.print('\t');
										manifest.println(kg.count_variants);
										}
										
									for(KeyGene kg: keyGenes) {
										Files.delete(kg.tmpVcfPath);
										}
									keyGenes.clear();
									
									if(ctx==null) break;
									prevCtg = ctx.getContig();
									}
						
					
								for(final GeneExtractorFactory.GeneExtractor ex: extractors)
									{
									final Map<GeneExtractorFactory.KeyAndGene,Set<String>> gene2values = ex.apply(ctx);
									
									if(gene2values.isEmpty()) continue;
									
									for(final GeneExtractorFactory.KeyAndGene keyAndGene :gene2values.keySet()) {
										final Set<String> values = gene2values.get(keyAndGene);
										if(values.isEmpty()) continue;
				
										KeyGene keyGene = keyGenes.stream().
													filter(KG->KG.extractor.equals(keyAndGene.getMethod())&& KG.key.equals(keyAndGene.getKey())).
													findFirst().orElse(null);
										if(keyGene==null) {
											keyGene = new KeyGene(keyAndGene.getMethod(),keyAndGene.getKey(),keyAndGene.getGene());
											
											keyGenes.add(keyGene);
											}
										final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
										vcb.rmAttribute(VepPredictionParser.getDefaultTag());
										vcb.rmAttribute(AnnPredictionParser.getDefaultTag());
										vcb.rmAttribute(BcfToolsPredictionParser.getDefaultTag());
										
											
										vcb.attribute(ex.getInfoTag(), new ArrayList<>(values));
				
										keyGene.write(header,vcb.make());
										
										
										final long num_files_opened = keyGenes.stream().filter(K->K.pw!=null).count();
										if(num_files_opened > this.max_open_files) {
											final KeyGene toClose = keyGenes.stream().
												filter(K->K.pw!=null).
												sorted((A,B)->Long.compare(A.lastModificationDate,B.lastModificationDate)).
												findFirst().
												orElse(null)
												;
											if(toClose!=null) {
												toClose.pw.flush();
												toClose.pw.close();
												toClose.pw= null;
												}
											}
										}
									}
						}// end of for(;;)
						
						}//end of in
			

					}// end of manifest
				}// end of archive
			return 0;
			}
		catch(final Throwable err) 
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			for(KeyGene kg:keyGenes) {
				if(kg.pw!=null) try {kg.pw.close();} catch(Throwable err) {}
				if(kg.tmpVcfPath!=null) try {Files.delete(kg.tmpVcfPath);} catch(Throwable err) {}
				}
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.list_extractors) {
			for(final String en: GeneExtractorFactory.getExtractorNames()) {
				System.out.println(en);
				}
			return 0;
			}
		
		
		try
			{
			return run(args);
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	 	
	
	public static void main(final String[] args)
		{
		new VcfGeneSplitter().instanceMainWithExit(args);
		}
	}
