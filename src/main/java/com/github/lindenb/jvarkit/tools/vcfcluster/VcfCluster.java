/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfcluster;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

/**
BEGIN_DOC

## Motivation

split CNV/SV vcf by clusters of size='x' for later parallelization (e.g: duphold, annotation, etc...)
All INFO/SVTYPE=BND are grouped in the same cluster using ID and INFO/MATEID

## Example

```bash
java -jar dist/jvarkit.jar -o

```

END_DOC
*/
@Program(name="vcfcluster",
	description="VCF",
	keywords={"vcf"},
	creationDate = "20260108",
	modificationDate = "20260108",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfCluster extends Launcher {
	private static Logger LOG=Logger.of(VcfCluster.class);
	@Parameter(names={"-o","--output"},description="output directory", required = true)
	private File outputDir = null;
	@Parameter(names={"-L","--length","-n"},description="max-length per cluster. " + DistanceParser.OPT_DESCRIPTION, converter = DistanceParser.StringConverter.class, splitter = NoSplitter.class)
	private int max_len = 1_000_000;
	@Parameter(names={"-seed","--seed"},description="random see: -1= current-time. Random is used to add a variant to a random available batch")
	private long random_seed = -1L;
	@Parameter(names={"-p","--prefix"},description="file prefix")
	private String prefix="split.";
	@Parameter(names={"-f","--force"},description="force writing")
	private boolean force_write = false;
	@Parameter(names={"--no-index"},description="do no write tbi index")
	private boolean disable_index = false;

	private static int ID_GENERATOR=0;
	
	/**
	 * Batch
	 */
	private  class Batch {
		final List<VariantContext> variants =new ArrayList<>();
		boolean bnd_flag = false;
	
		long getLengthOnReference() {
			return VcfCluster.getLengthOnReference(variants);
			}
		
		void save(final VariantContextWriterBuilder vcw,final VCFHeader header) throws IOException {
			// special, save BND in the same cluster
			if(bnd_flag && getLengthOnReference()> VcfCluster.this.max_len) {
				Batch subbatch = null;
				while(!variants.isEmpty()) {
					final VariantContext ctx = variants.remove(0);
					if(subbatch==null) subbatch = new Batch();
					subbatch.bnd_flag=true;
					subbatch.variants.add(ctx);
					int i=0;
					// loop , always, check that BND linked together are present in the same cluster
					while(i< this.variants.size()) {
						final VariantContext ctx2 = this.variants.get(i);
						final String ctx2_id = ctx2.hasID()?ctx2.getID():null;
						final String ctx2_mate_id = ctx2.hasAttribute("MATEID")?ctx2.getAttributeAsString("MATEID", ""):null;
						
						for(VariantContext ctx3 : subbatch.variants) {
							final String ctx3_id = ctx3.hasID()?ctx3.getID():null;
							final String ctx3_mate_id = ctx3.hasAttribute("MATEID")?ctx3.getAttributeAsString("MATEID", ""):null;
							if(
								(!StringUtils.isBlank(ctx2_id) && ctx2_id.equals(ctx3_mate_id)) ||
								(!StringUtils.isBlank(ctx3_id) && ctx3_id.equals(ctx2_mate_id))
								) {
								 subbatch.variants.add(ctx2);
								 this.variants.remove(i);
								}
							else
								{
								i++;
								}
							}
						}
					if(subbatch!=null && subbatch.getLengthOnReference()>=VcfCluster.this.max_len) {
						save0(vcw,header,subbatch.variants);
						subbatch = null;
						}
					}
				if(subbatch!=null) {
					save0(vcw,header,subbatch.variants);
					}
				}
			else
				{
				save0(vcw,header,this.variants);
				}
			}
		}

	private static long getLengthOnReference(List<VariantContext> variants) {
		return variants.stream().mapToLong(V->V.getLengthOnReference()).sum();
		}
	
	private void save0(final VariantContextWriterBuilder vcw,final VCFHeader header,final List<VariantContext> variants) throws IOException {
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
		this.outputDir.mkdirs();
		final File filename = new File(outputDir,String.format("%s%d.L%s.N%d.vcf.gz", 
				this.prefix,
				(++ID_GENERATOR),
				String.valueOf(getLengthOnReference(variants)),
				variants.size()
				));
		if(filename.exists() && !force_write) {
			throw new IOException("File exists : "+filename);
			}
		LOG.warn("Save "+filename);
		
		Collections.sort(variants, new VariantContextComparator(dict));
		vcw.setOutputFile(filename);
		vcw.setOutputFileType(VariantContextWriterBuilder.determineOutputTypeFromFile(filename));
		try( VariantContextWriter  w=vcw.build()) {
			w.writeHeader(header);
			for(VariantContext vc: variants) {
				w.add(vc);
				}
			}
		}
	
@Override
public int doWork(final List<String> args) {
	try {
		final Random rand = new Random(random_seed==-1L?System.currentTimeMillis():random_seed);
		final String input =  super.oneFileOrNull(args);
		try(VCFIterator iter = input==null?new VCFIteratorBuilder().open(stdin()):new VCFIteratorBuilder().open(input)) {
			final VCFHeader header = iter.getHeader();
			JVarkitVersion.getInstance().addMetaData(this, header);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final VariantContextWriterBuilder vcw= new VariantContextWriterBuilder();
			vcw.setReferenceDictionary(dict);
			if(!disable_index) {
				final TabixIndexCreator tbi=new TabixIndexCreator(dict, TabixFormat.VCF);
				vcw.setIndexCreator(tbi);
				vcw.setOption(Options.INDEX_ON_THE_FLY);
				}
			vcw.unsetBuffering();
			
			final List<Batch> batches = new ArrayList<>();
			Batch bnd_batch=null;
			for(;;) {
				final VariantContext ctx = iter.hasNext()?iter.next():null;
				
				if(ctx==null) {
					for(Batch b:batches) {
						b.save(vcw,header);
						}
					if(bnd_batch!=null) {
						bnd_batch.save(vcw,header);
						bnd_batch = null;
						}
					break;
					}
				
				if(ctx.hasAttribute(VCFConstants.SVTYPE) && 
					ctx.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")) {
					if(bnd_batch==null) {
						bnd_batch = new Batch();
						bnd_batch.bnd_flag = true;
						}
					
					bnd_batch.variants.add(ctx);
					continue;
					}
				
				final int len = ctx.getLengthOnReference();

				if(len >= max_len) {
					final Batch b=new Batch();
					b.variants.add(ctx);
					b.save(vcw,header);
					continue;
					}
				
			
				final List<Batch> copy=new ArrayList<>(batches);
				copy.removeIf(B->B.getLengthOnReference()+len > max_len);
				if(!copy.isEmpty()) {
					Collections.shuffle(copy,rand);
					copy.get(0).variants.add(ctx);
					}
				else
					{
					final Batch b=new Batch();
					b.variants.add(ctx);
					batches.add(b);
					}
				int i=0;
				while(i< batches.size()) {
					if( (double)batches.get(i).getLengthOnReference() >= max_len*0.9 || batches.get(i).variants.size()>100) {
						batches.get(i).save(vcw,header);
						batches.remove(i);
						}
					else
						{
						i++;
						}
					}
				}
			}
	   
		return 0;
		}
	catch(final Throwable err) {
		err.printStackTrace();
		return -1;
		}
	}

public static void main(final String[] args) {
	new VcfCluster().instanceMain(args);
	}
}

