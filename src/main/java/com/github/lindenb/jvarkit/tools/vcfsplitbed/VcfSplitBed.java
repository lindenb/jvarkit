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
import java.io.OutputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.MultiIntervalVariantIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

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
	@Parameter(names={"--prefix"},description="prefix each output VCF file with this string")
	private String prefix="";
	@Parameter(names={"--disable-hash-directory","--dhd"},description="disable default which is to save each file in a checksum-based directory-a-la-nextflow to avoid a large number of files in the same directory.")
	private boolean disable_hash_dir = false;

	
	@Override
	public int doWork(final List<String> args) {
    	
		try
			{
			final Map<String,List<Locatable>> id2intervals = new HashMap<>();
			final Path vcfPath = Paths.get(oneAndOnlyOneFile(args));
			try(VCFFileReader vcfFileReader = new VCFFileReader(vcfPath,true)) {
				final VCFHeader header0  = vcfFileReader.getHeader();
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
						final Locatable rgn = new SimpleInterval(rec.getContig(), rec.getStart(), rec.getEnd());
						List<Locatable> gene_set = id2intervals.get(intervalName);
						if(gene_set==null) {
							gene_set = new ArrayList<>();
							id2intervals.put(intervalName, gene_set);
							}
						gene_set.add(rgn);
						}
					} // end loop bed
					LOG.info("initial number of identifiers  =  "+id2intervals.size());
					LOG.info("initial number of 'intervals'  =  "+id2intervals.values().stream().mapToLong(G->G.stream().count()).sum());
					
					
					try(ArchiveFactory archiveFactory = ArchiveFactory.open(this.outputFile)) {
						for(String id: id2intervals.keySet()) {
							try(CloseableIterator<VariantContext> r = MultiIntervalVariantIterator.query(vcfFileReader, id2intervals.get(id))) {
								if(!r.hasNext()) continue;
								final String md5 = StringUtils.md5(id);
								final String parentDir = md5.substring(0,2) + File.separatorChar + md5.substring(2);
								final String filename0 =
										(this.disable_hash_dir?"":parentDir + File.separator)+
										this.prefix+
										id.replaceAll("[/\\:_]+", "_") + ".vcf.gz";
								
								final VCFHeader h2 = new VCFHeader(header0);
								h2.addMetaDataLine(new VCFHeaderLine("VcfSplitBed.Name",String.valueOf(id)));
								JVarkitVersion.getInstance().addMetaData(VcfSplitBed.this, h2);
								try(OutputStream os = archiveFactory.openOuputStream(filename0))  {
									try(VariantContextWriter vcw = VCFUtils.createVariantContextWriterToOutputStream(os)) {
										vcw.writeHeader(h2);
										while(r.hasNext()) {
											vcw.add(r.next());
											}
										} // vcw
									} // output
							} // end iterator
						} // end loop geneid
					}//end archive
			
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	 	
	
	public static void main(final String[] args)
		{
		new VcfSplitBed().instanceMainWithExit(args);
		}
	}
