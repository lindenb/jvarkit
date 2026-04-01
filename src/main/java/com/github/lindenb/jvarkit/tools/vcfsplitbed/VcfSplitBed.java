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
package com.github.lindenb.jvarkit.tools.vcfsplitbed;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC

### Example

```
$ java -jar dist/jvarkit.jar vcfsplitbed --bed input.ned src/test/resources/rotavirus_rf.vcf.gz 
```

# see also

 * vcfgenesplitter

END_DOC
*/
@Program(
		name="vcfsplitbed",
		description="Split VCF by BED/Region Name",
		creationDate="20190619",
		modificationDate="20260401",
		jvarkit_amalgamion = true,
		keywords= {"vcf","sliding","window"}
		)
public class VcfSplitBed
	extends Launcher
	{
	private static final Logger LOG = Logger.of(VcfSplitBed.class);
		
	@Parameter(names={"-o","--output"},description= ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-B","--bed"},description="User bed fourth column is used as the key to aggregate the variants. Empty is replace with {chrom}_{start}_{end} ",required = true)
	private Path userBed=null;
	@Parameter(names={"--open-max"},description="Maximum number of opened VCF writers at the same time.")
	private int max_open_files = 100;
	@Parameter(names={"-n","--min-variant"},description="Minimum number of variants required to write a vcf. don't write if num(variant) < 'x' ")
	private int min_number_of_ctx = 1;
	@Parameter(names={"-M","--max-variant"},description="Maximum number of variants required to write a vcf. don't write if num(variant) > 'x' . '<=0' is ignore")
	private int max_number_of_ctx = -1;
	@Parameter(names={"--prefix"},description="prefix each output VCF file with this string")
	private String prefix="";
	@Parameter(names={"--disable-hash-directory","--dhd"},description="disable default which is to save each file in a checksum-based directory-a-la-nextflow to avoid a large number of files in the same directory.")
	private boolean disable_hash_dir = false;



	private class IntervalDefinition {

		final String intervalName;
		int count_variants = 0;
		Path tmpVcfPath = null;
		PrintWriter pw = null;
		VCFEncoder vcfEncoder = null;
		long lastModificationDate = 0L;
		

		
		IntervalDefinition(final String intervalName)  throws IOException {
			this.intervalName = intervalName;
			}
		@Override
		public int hashCode() {
			return this.intervalName.hashCode();
			}
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof IntervalDefinition)) return false;
			final IntervalDefinition kg = IntervalDefinition.class.cast(obj);
			return this.intervalName.equals(kg.intervalName);
			}
	
		public void write(final VCFHeader header,final VariantContext ctx) throws IOException {
			if(this.pw!=null) {
				//nothing
				}
			else if(this.tmpVcfPath ==null ) {
				LOG.info("Opening VCF for "+this.intervalName);
				this.tmpVcfPath =  Files.createTempFile("tmp.", ".vcf");
				this.pw = new PrintWriter(Files.newBufferedWriter(this.tmpVcfPath, StandardOpenOption.APPEND));
				final VCFHeader h2 = new VCFHeader(header);
				h2.addMetaDataLine(new VCFHeaderLine("VcfSplitBed.Name",String.valueOf(this.intervalName)));
				JVarkitVersion.getInstance().addMetaData(VcfSplitBed.this, h2);
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
			this.vcfEncoder.write(this.pw,ctx);
			this.pw.println();
			}
		@Override
		public String toString() {
			return  this.intervalName;
			}
		void mayClose() {
			if(this.pw!=null) {
				try {
				this.pw.flush();
				this.pw.close();
				} catch(Throwable err) {
					LOG.error(err);
					}
				this.pw=null;
				}
			}
		}
	
	
		
	public VcfSplitBed()
		{
		}
	
	
	
	
	@Override
	public int doWork(final List<String> args) {
		
    	
    	if(this.min_number_of_ctx<=0) {
    		LOG.error("Bad minimum number of variants");
    		return -1;
    	}
    	
    	
		final Map<String,IntervalDefinition> gene2keygene = new HashMap<>();
    	
		try
			{
			final IntervalTreeMap<Set<String>> interval2geneid = new IntervalTreeMap<>();
			try(VCFIterator r=super.openVCFIterator(oneFileOrNull(args))) {
				final VCFHeader header0  = r.getHeader();
				/* load bed */
				final SAMSequenceDictionary dict = header0.getSequenceDictionary();
				try(BedLineReader br= new BedLineReader(this.userBed)) {
					if(dict!=null) br.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					while(br.hasNext()) {
						final BedLine rec = br.next();
						if(rec==null) continue;
						String intervalName = rec.getOrDefault(3, "");
						if(StringUtils.isBlank(intervalName)) {
							intervalName = rec.getContig()+"_"+rec.getStart()+"_"+rec.getEnd();
							}
						final Interval rgn = new Interval(rec.getContig(), rec.getStart(), rec.getEnd());
						Set<String> gene_set = interval2geneid.get(rgn);
						if(gene_set==null) {
							gene_set = new HashSet<>();
							interval2geneid.put(rgn, gene_set);
							}
						gene_set.add(intervalName);
						}
					} // end loop bed
					LOG.info("initial number of intervals  =  "+interval2geneid.size());
					LOG.info("initial number of 'aggregates'  =  "+interval2geneid.values().stream().flatMap(G->G.stream()).collect(Collectors.toSet()).size());
					
					
					while(r.hasNext()) {
						final VariantContext ctx = r.next();
						for(final String gene_id: interval2geneid.getOverlapping(ctx).stream().flatMap(C->C.stream()).collect(Collectors.toSet())) {
							IntervalDefinition keyGene = gene2keygene.get(gene_id);
							if(keyGene==null) {
								keyGene = new IntervalDefinition(gene_id);
								gene2keygene.put(gene_id,keyGene);
								}
							
							if ( this.max_number_of_ctx!=-1 && keyGene.count_variants > this.max_number_of_ctx ) {
								continue;
								}
							
							keyGene.write(header0,ctx);
							
							final long num_files_opened = gene2keygene.values().stream().filter(K->K.pw!=null).count();
							if(num_files_opened > this.max_open_files) {
								final IntervalDefinition toClose = gene2keygene.values().stream().
									filter(K->K.pw!=null).
									sorted((A,B)->Long.compare(A.lastModificationDate,B.lastModificationDate)).
									findFirst().
									orElse(null)
									;
								toClose.mayClose();
								
								} // end if num_file
							}// end loop over gene_id
						}// end while
					
				
				}//end vcf iterator
			
			try(ArchiveFactory archiveFactory = ArchiveFactory.open(this.outputFile)) {
				for(IntervalDefinition kg: gene2keygene.values()) {
					kg.mayClose();
					
					if ( kg.count_variants < this.min_number_of_ctx )  {
						LOG.info("skipping "+kg+" because there are not enough variants. N="+kg.count_variants+"<"+this.min_number_of_ctx);
						continue;
						}
					if ( this.max_number_of_ctx!=-1 && kg.count_variants > this.max_number_of_ctx ) {
						LOG.info("skipping "+kg+" because there are too many variants. N="+kg.count_variants+">"+this.max_number_of_ctx);
						continue;
						}
					
					final String md5 = StringUtils.md5(kg.intervalName);
					final String parentDir = md5.substring(0,2) + File.separatorChar + md5.substring(2);
					final String filename0 =
							(this.disable_hash_dir?"":parentDir + File.separator)+
							this.prefix+
							kg.intervalName.replaceAll("[/\\:_]+", "_") + ".vcf.gz";
					
					
					try(final BlockCompressedOutputStream os = new BlockCompressedOutputStream(archiveFactory.openOuputStream(filename0),(Path)null)) {
						IOUtils.copyTo(kg.tmpVcfPath, os);
						os.flush();
						}
					Files.delete(kg.tmpVcfPath);
					} //end loop over gene
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			for(IntervalDefinition kg: gene2keygene.values()) {
				try {
					kg.mayClose();
					Files.delete(kg.tmpVcfPath);
					}
				catch(IOException err) {
					
					}
				}
			}
		}
	 	
	
	public static void main(final String[] args)
		{
		new VcfSplitBed().instanceMainWithExit(args);
		}
	}
