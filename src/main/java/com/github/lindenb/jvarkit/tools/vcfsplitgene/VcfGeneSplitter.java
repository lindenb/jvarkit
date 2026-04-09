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
package com.github.lindenb.jvarkit.tools.vcfsplitgene;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
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
#chrom	POS	key	path	Count_Variants
RF01	969		ANN/GeneId	Gene_18_3284	2c/8fb9d2539e3f30d1d9b06f9ec54c4c/Gene_18_3284.vcf.gz	1
RF02	250		ANN/GeneId	Gene_1621_1636	4e/4897c51fe2dd067a8b75c19f111477/Gene_1621_1636.vcf.gz	5
RF02	250		ANN/GeneId	UniProtKB/Swiss-Prot:P12472	74/ca4273c3d5803c5865891c808234da/UniProtKB_Swiss-Prot:P12472.vcf.gz	5
RF03	1220		ANN/GeneId	Gene_50_2557	23/6b59cfe4fdd33a5f4feeb55521dd34/Gene_50_2557.vcf.gz	8
RF04	886		ANN/GeneId	Gene_9_2339	b3/4bda8d8502e64e442fce077e45ded6/Gene_9_2339.vcf.gz	7
RF05	40		ANN/GeneId	Gene_32_1507	b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz	6
RF06	516		ANN/GeneId	Gene_23_1216	6f/8472e9f192c92bf46e4893b2367b7e/Gene_23_1216.vcf.gz	5
RF07	97		ANN/GeneId	Gene_0_1073	3c/513d82eaea18447dd5f621f92b40e6/Gene_0_1073.vcf.gz	4
RF08	925		ANN/GeneId	Gene_0_1058	84/977eac8cdef861cbd3109209675d21/Gene_0_1058.vcf.gz	2
RF09	293		ANN/GeneId	Gene_0_1061	db/aee9cc8f5c9c3d39c7af4cec63b7a5/Gene_0_1061.vcf.gz	3
RF10	45		ANN/GeneId	Gene_41_568	b0/133c483f0ea676f8d29ab1f2daee5d/Gene_41_568.vcf.gz	3
RF11	73		ANN/GeneId	Gene_20_616	59/0fd5c1e8d6d60a986a0021fe357514/Gene_20_616.vcf.gz	1
RF11	73		ANN/GeneId	Gene_78_374	83/bc905cf311428ab80ce59aaf503838/Gene_78_374.vcf.gz	1


```

## See also

* vcfwindowsplitter

END_DOC
*/
@Program(
		name="vcfgenesplitter",
		description="Split VCF+VEP by gene/transcript.",
		creationDate = "20160310",
		modificationDate="20260409",
		keywords= {"genes","vcf"},
		jvarkit_amalgamion =  true,
		menu="VCF Manipulation"
		)
public class VcfGeneSplitter
	extends Launcher
	{
	private static final Logger LOG = Logger.of(VcfGeneSplitter.class);
	

    private static class GeneName
	    {
	    final String gene_id;
	    final String label;
	    final String extractorName;
	    GeneName(final String gene_id,final String label,final String type)
	            {
	            this.gene_id=gene_id;
	            this.label=StringUtils.isBlank(label)?".":label;
	            this.extractorName=type;
	            }
	    @Override
	    public int hashCode()
	            {
	            final int prime = 31;
	            int result = 1;
	            result = prime * result +  gene_id.hashCode();
	            result = prime * result +  extractorName.hashCode();
	            return result;
	            }
	    @Override
	    public boolean equals(final Object o)
	            {
	            if (this == o) return true;
	            if (o == null) return false;
	            if (getClass() != o.getClass()) return false;
	            final GeneName g=(GeneName)o;
	            return gene_id.equals(g.gene_id) && extractorName.equals(g.extractorName);
	            }
	    @Override
	    public String toString() {
	            return  gene_id+"("+extractorName+")";
	            }
	
	    }

	
    private class Call implements Comparable<Call>
	    {
    	GeneName gene;
	    VariantContext ctx;
	
	    Call() {
	    	this(null,null);
	    	}
	    Call(GeneName gene,VariantContext ctx) {
	    	this.gene = gene;
	    	this.ctx = ctx;
	    	}
	
	    String getContig()
	            {
	            return ctx.getContig();
	            }
	
	    @Override
	    public int compareTo(final Call o) {
	            int i=  this.getContig().compareTo(o.getContig());
	            if(i!=0) return i;
	            i= this.gene.gene_id.compareTo(o.gene.gene_id);
	            if(i!=0) return i;
	            i= this.gene.extractorName.compareTo(o.gene.extractorName);
	            return i;
	            }
	
	    public int compare2(final Call C2) {
	            int i= this.compareTo(C2);
	            if(i!=0) return i;
	            i =  this.ctx.getContig().compareTo(C2.ctx.getContig());
	            if(i!=0) return i;
	            i =  Integer.compare(this.ctx.getStart(),C2.ctx.getStart());
	            if(i!=0) return i;
	            i =  this.ctx.getReference().compareTo(C2.ctx.getReference());
	            return i;
	            }
	
	
	    }
    private class CallCodec 
    extends AbstractDataCodec<Call>
	    {
	    final VCFHeader header;
	    private final VCFCodec vCodec;
	    private final VCFEncoder vcfEncoder;
	
	    CallCodec(final VCFHeader header) {
	            this.header= header;
	            this.vCodec = new VCFCodec();
	            this.vcfEncoder = new VCFEncoder(header, false, false);
	            this.vCodec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
	            }
	    @Override
	    public void encode(final DataOutputStream dos,final Call c)
	                    throws IOException
	            {
	            dos.writeUTF(c.gene.gene_id);
	            dos.writeUTF(c.gene.label);
	            dos.writeUTF(c.gene.extractorName);
	            writeString(dos, this.vcfEncoder.encode(c.ctx));
	            }
	
	    @Override
	    public Call decode(final DataInputStream dis) throws IOException
	            {
	            final String gene_id;
	            try {
	            	gene_id=dis.readUTF();
	            } catch (final EOFException e) {
	                    return null;
	                    }
	            final String label=dis.readUTF();
	            final String extractor=dis.readUTF();
	            final Call c= new Call();
	            c.gene=new GeneName(gene_id,label, extractor);
	            c.ctx = this.vCodec.decode(readString(dis));
	            return c;
	            }
	    @Override
	    public CallCodec clone() {
	            return new CallCodec(this.header);
	            }
	    	}
	
	@Parameter(names={"-o","--output"},description= ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-m","--manifest"},description="Manifest BED file output containing chrom/POS of each gene")
	private Path manifestFile = null;
	@Parameter(names={"-l","--list"},description= "list all available extractors", help=true)
	private boolean list_extractors = false;
	@Parameter(names={"-e","-E","--extractors"},description=GeneExtractorFactory.OPT_DESC)
	private String extractorsNames="ANN/GeneId VEP/GeneId";
	@Parameter(names={"--ignore-filtered"},description="Ignore FILTERED variant")
	private boolean ignoreFiltered = false;
	@Parameter(names={"--prefix"},description="prefix each output VCF file with this string")
	private String prefix="";
	@Parameter(names={"--disable-hash-directory","--dhd"},description="disable default which is to save each file in a checksum-based directory-a-la-nextflow to avoid a large number of files in the same directory.")
	private boolean disable_hash_dir = false;
    @ParametersDelegate
    private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

	
	
	public VcfGeneSplitter()
		{
		
		}

	
	private int run(final List<String> args) {
		SortingCollection<Call> sortingCollection=null;
		final VCFHeader vcfHeader ;
		try {
		try(VCFIterator iterator = super.openVCFIterator(oneFileOrNull(args))) {
			vcfHeader = iterator.getHeader();
			final GeneExtractorFactory geneExtractorFactory = new GeneExtractorFactory(vcfHeader);
			final List<GeneExtractorFactory.GeneExtractor> extractors = geneExtractorFactory.parse(this.extractorsNames);
			if(extractors.isEmpty()) {
				LOG.warn("No extractor defined!");
				return -1;
				}
			
            sortingCollection =SortingCollection.newInstance(
                    Call.class,
                    new CallCodec(vcfHeader),
                    (C1,C2)->C1.compare2(C2),
                    this.writingSortingCollection.getMaxRecordsInRam(),
                    this.writingSortingCollection.getTmpPaths()
                    );
            sortingCollection.setDestructiveIteration(true);
            while(iterator.hasNext()) {
            	final VariantContext ctx = iterator.next();
				
				if(this.ignoreFiltered && ctx!=null && ctx.isFiltered()) continue;
				for(final GeneExtractorFactory.GeneExtractor ex: extractors)
					{
					final Map<GeneExtractorFactory.KeyAndGene,Set<String>> gene2values = ex.apply(ctx);
					for(final GeneExtractorFactory.KeyAndGene keyAndGene :gene2values.keySet()) {
						final Set<String> values = gene2values.get(keyAndGene);
						final GeneName geneName = new GeneName(keyAndGene.getKey(), StringUtils.ifBlank(keyAndGene.getGene(),"."), keyAndGene.getMethod());
						
						
						final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
						vcb.rmAttribute(VepPredictionParser.getDefaultTag());
						vcb.rmAttribute(AnnPredictionParser.getDefaultTag());
						vcb.rmAttribute(BcfToolsPredictionParser.getDefaultTag());
						
						if(ex.hasInfoTag() && !values.isEmpty()) {
							vcb.attribute(ex.getInfoTag(), new ArrayList<>(values));
							}
						sortingCollection.add(new Call(geneName,vcb.make()));
						} // end genes
	            	}//end extractors
            	}// end iterator VCF
			}
			sortingCollection.doneAdding();
			try(CloseableIterator<Call> iter = sortingCollection.iterator()) {
				try(ArchiveFactory archiveFactory = ArchiveFactory.open(this.outputFile)) {
					try(PrintWriter manifest = new PrintWriter(this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(manifestFile))) {
						manifest.println("#chrom\tstart\tend\tsplitter\tgene\tkey\tpath\tCount_Variants");
						
						
						GeneName prevGeneName =   null;
						VariantContextWriter currentWriter = null;
						int min_pos=-1;
						int max_pos=-1;
						int count_variants=0;
						String currContig=null;
						String filename0 = null;
						for(;;)
							{
							final Call curr = iter.hasNext()?iter.next():null;
							
							if(curr==null || !curr.gene.equals(prevGeneName)) {
								if(currentWriter!=null) {
									manifest.print(currContig);
									manifest.print('\t');
									manifest.print(min_pos-1);
									manifest.print('\t');
									manifest.print(max_pos);
									manifest.print('\t');
									manifest.print(prevGeneName.extractorName);
									manifest.print('\t');
									manifest.print(prevGeneName.label);
									manifest.print('\t');
									manifest.print(prevGeneName.gene_id);
									manifest.print('\t');
									manifest.print(
										archiveFactory.isTarOrZipArchive()?
										filename0:
										this.outputFile.resolve(filename0).toAbsolutePath().toString()
										);
									manifest.print('\t');
									manifest.println(count_variants);
								
									currentWriter.close();
									currentWriter=null;
									
									if(curr==null) break;
									
								
									}
								min_pos = curr.ctx.getStart();
								max_pos= curr.ctx.getEnd();
								count_variants = 0;
								filename0 = null;
								
								prevGeneName = curr.gene;
								
								
								
								final String md5 = StringUtils.md5(curr.getContig()+":"+curr.gene.extractorName+":"+curr.gene.gene_id);
								final String parentDir = md5.substring(0,2) + File.separatorChar + md5.substring(2);
								filename0 =
										(this.disable_hash_dir?"":parentDir + File.separator)+
										this.prefix+
										(curr.getContig()+"_"+curr.gene.extractorName+"_"+curr.gene.label+"_"+curr.gene.gene_id).replaceAll("[/\\:_]+", "_") + ".vcf.gz";
								final VCFHeader h2 = new VCFHeader(vcfHeader);
								h2.addMetaDataLine(new VCFHeaderLine("VcfGeneSplitter.GeneId",String.valueOf(curr.gene.gene_id)));
								h2.addMetaDataLine(new VCFHeaderLine("VcfGeneSplitter.GeneName",String.valueOf(curr.gene.label)));
								h2.addMetaDataLine(new VCFHeaderLine("VcfGeneSplitter.Extractor",String.valueOf(curr.gene.extractorName)));
								JVarkitVersion.getInstance().addMetaData(VcfGeneSplitter.this, h2);
								
								final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
								vcwb.setCreateMD5(false);
								vcwb.setReferenceDictionary(vcfHeader.getSequenceDictionary());
								vcwb.clearOptions();
								vcwb.setOutputStream(new BlockCompressedOutputStream(archiveFactory.openOuputStream(filename0),(Path)null));
								currentWriter =vcwb.build();
								currentWriter.writeHeader(h2);
								}
							currentWriter.add(curr.ctx);
							currContig = curr.ctx.getContig();
							min_pos = Math.min(min_pos,curr.ctx.getStart());
							max_pos= Math.max(max_pos,curr.ctx.getEnd());
							count_variants++;
							}
						manifest.flush();
						}//end manifest
					}
				}
			sortingCollection.cleanup();
			sortingCollection=null;
			return 0;
			}
		catch(final Throwable err) 
			{
			LOG.error(err);
			return -1;
			}
		finally {
			if(sortingCollection!=null) try {
				sortingCollection.cleanup();
			} catch(Throwable err) {
				
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
