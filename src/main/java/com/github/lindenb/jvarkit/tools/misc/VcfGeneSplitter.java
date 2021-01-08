/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

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
		modificationDate="20190614",
		keywords= {"genes","vcf"}
		)
public class VcfGeneSplitter
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfGeneSplitter.class).make();
	
	private static class KeyGene
		implements Comparable<KeyGene>
		{
		final String key;
		final String gene;
		KeyGene(final String key,final String gene) {
			this.key = key;
			this.gene = StringUtils.isBlank(gene)?".":gene;
			}
		@Override
		public int hashCode() {
			return key.hashCode();
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof KeyGene)) return false;
			return this.compareTo(KeyGene.class.cast(obj))==0;
			}
		@Override
		public int compareTo(final KeyGene o) {
			return this.key.compareTo(o.key);
			}
		@Override
		public String toString() {
			return this.key;
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


	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	private static class KeyAndLine {
		final KeyGene keyAndGene;
		final String splitter;
		final String ctx;
		KeyAndLine(final KeyGene keyAndGene,final String splitter,final String ctx) {
			this.keyAndGene = keyAndGene;
			this.splitter = splitter;
			this.ctx = ctx;
		}
	}
	
	private static class KeyAndLineComparator1
		implements Comparator<KeyAndLine>
		{
		@Override
		public int compare(final KeyAndLine o1, final KeyAndLine o2) {
			int i = o1.keyAndGene.compareTo(o2.keyAndGene);
			if(i!=0) return i;
			i = o1.splitter.compareTo(o2.splitter);
			return i;
			}
		}
	private static class KeyAndLineComparator2
		extends KeyAndLineComparator1
		{
		final CharSplitter tab = CharSplitter.TAB;
		@Override
		public int compare(final KeyAndLine o1, final KeyAndLine o2) {
			int i = super.compare(o1, o2);
			if(i!=0) return i;
			final String tokens1[] = this.tab.split(o1.ctx,5);
			final String tokens2[] = this.tab.split(o2.ctx,5);
			if(!tokens1[0].equals(tokens2[0])) {
				throw new IllegalStateException("not same contig???");
			}
			i = Integer.parseInt(tokens1[1]) - Integer.parseInt(tokens2[1]);
			if(i!=0) return i;
			i = Allele.create(tokens1[3],true).compareTo( Allele.create(tokens2[3],true));
			return i;
			}
		}
	
	private static class KeyAndLineCodec extends AbstractDataCodec<KeyAndLine>
		{
		@Override
		public KeyAndLine decode(final DataInputStream dis) throws IOException {
			String k;
			try {
				k=dis.readUTF();
			} catch(EOFException err) { return null;}
			final String g = dis.readUTF();
			final String splt = dis.readUTF();
			final String v = AbstractDataCodec.readString(dis);
			return new KeyAndLine(new KeyGene(k, g),splt, v);
		}
		@Override
		public void encode(final DataOutputStream dos,final  KeyAndLine object) throws IOException {
			dos.writeUTF(object.keyAndGene.key);
			dos.writeUTF(object.keyAndGene.gene);
			dos.writeUTF(object.splitter);
			AbstractDataCodec.writeString(dos, object.ctx);
			}
		@Override
		public AbstractDataCodec<KeyAndLine> clone() {
			return new KeyAndLineCodec();
			}
		}

		
	public VcfGeneSplitter()
		{
		
		}
	
	
	private int run(final List<String> args) {
		SortingCollection<KeyAndLine> sortingcollection=null;
		BufferedReader in = null;
		FileOutputStream fos = null;
		CloseableIterator<KeyAndLine> iter=null;
		ArchiveFactory archiveFactory = null;
		PrintWriter manifest = null;
		try {
			final Path tmpVcf = Files.createTempFile("tmp.", ".vcf.gz");
			
			archiveFactory = ArchiveFactory.open(this.outputFile);
			
			manifest = new PrintWriter(this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(manifestFile));
			manifest.println("#chrom\tstart\tend\tsplitter\tgene\tkey\tpath\tCount_Variants");

			in = super.openBufferedReader(oneFileOrNull(args));
			final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(in);
			
			
			
			final GeneExtractorFactory geneExtractorFactory = new GeneExtractorFactory(cah.header);
			
			final List<GeneExtractorFactory.GeneExtractor> extractors = geneExtractorFactory.parse(this.extractorsNames);
			if(extractors.isEmpty()) {
				LOG.warn("No extractor defined!");
				return -1;
				}
			
			
			final VCFEncoder vcfEncoder = new VCFEncoder(cah.header, false, false);
			
			
			// read variants
			final ProgressFactory.Watcher<VariantContext> progess= ProgressFactory.newInstance().dictionary(cah.header).logger(LOG).build();
			String prevCtg = null;
			for(;;)
				{
				final String line = in.readLine();
				
				final VariantContext ctx = (line==null?null:cah.codec.decode(line));
				if(ctx!=null) progess.apply(ctx);
				if(this.ignoreFiltered && ctx.isFiltered()) continue;
				
				if(ctx==null || !ctx.getContig().equals(prevCtg))
					{
					if(sortingcollection!=null) {
						sortingcollection.doneAdding();
						iter = sortingcollection.iterator();
						final EqualRangeIterator<KeyAndLine> eqiter = new EqualRangeIterator<>(iter,
							(o1,o2)->o1.keyAndGene.compareTo(o2.keyAndGene)
							);
						
						while(eqiter.hasNext()) {
							final List<KeyAndLine> buffer = eqiter.next();
							
							if ( buffer.size() < this.min_number_of_ctx ) continue;
							if ( this.max_number_of_ctx > 0 && buffer.size() > this.max_number_of_ctx ) continue;

							
							final KeyAndLine first = buffer.get(0);
							
							final VariantContextWriter out = VCFUtils.createVariantContextWriterToPath(tmpVcf);
							final VCFHeader header2=addMetaData(new VCFHeader(cah.header));
							header2.addMetaDataLine(new VCFHeaderLine("GtfFileSplitter.Name",String.valueOf(first.keyAndGene.key)));
							header2.addMetaDataLine(new VCFHeaderLine("GtfFileSplitter.Gene",String.valueOf(first.keyAndGene.gene)));
							out.writeHeader(header2);
							int minPos=Integer.MAX_VALUE;
							int maxPos=0;
							for(final KeyAndLine kl:buffer) {
								final VariantContext ctx2 = cah.codec.decode(kl.ctx);
								minPos = Math.min(ctx2.getStart(), minPos);
								maxPos = Math.max(ctx2.getEnd(), maxPos);
								out.add(ctx2);
								}
							out.close();
							
							final String md5 = StringUtils.md5(prevCtg+":"+first.splitter+":"+first.keyAndGene.key);
							final String filename =  md5.substring(0,2) + File.separatorChar + md5.substring(2) + File.separator+buffer.get(0).keyAndGene.key.replaceAll("[/\\:]", "_") + ".vcf.gz";
							
							
							final OutputStream os = archiveFactory.openOuputStream(filename);
							IOUtils.copyTo(tmpVcf, os);
							os.flush();
							os.close();
							
							manifest.print(prevCtg);
							manifest.print('\t');
							manifest.print(minPos-1);
							manifest.print('\t');
							manifest.print(maxPos);
							manifest.print('\t');
							manifest.print(first.splitter);
							manifest.print('\t');
							manifest.print(first.keyAndGene.gene);
							manifest.print('\t');
							manifest.print(first.keyAndGene.key);
							manifest.print('\t');
							manifest.print((archiveFactory.isTarOrZipArchive()?"":this.outputFile.toString()+File.separator)+filename);
							manifest.print('\t');
							manifest.println(buffer.size());
							}
						
						eqiter.close();
						iter.close();
						
					
						}
					sortingcollection = null;
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
						
						final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
						vcb.rmAttribute(VepPredictionParser.getDefaultTag());
						vcb.rmAttribute(AnnPredictionParser.getDefaultTag());
						
							
						vcb.attribute(ex.getInfoTag(), new ArrayList<>(values));
						
						
						
						if(sortingcollection==null) {
							sortingcollection = SortingCollection.newInstance(
									KeyAndLine.class,
									new KeyAndLineCodec(),
									new KeyAndLineComparator2(),
									this.writingSortingCollection.getMaxRecordsInRam(),
									this.writingSortingCollection.getTmpPaths()
									);
							sortingcollection.setDestructiveIteration(true);
							}
						
						sortingcollection.add(
							new KeyAndLine(
									new KeyGene(keyAndGene.getKey(),keyAndGene.getGene()),
									ex.getName(), 
									vcfEncoder.encode(vcb.make())
							));
						}
					}
				}
			progess.close();
			manifest.flush();
			manifest.close();
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
			CloserUtil.close(iter);
			if(sortingcollection!=null) sortingcollection.cleanup();
			CloserUtil.close(in);
			CloserUtil.close(fos);
			CloserUtil.close(manifest);
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
