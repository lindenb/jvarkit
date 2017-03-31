package com.github.lindenb.jvarkit.tools.vcfcomposite;

import java.io.File;
import java.io.IOException;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="vcfcomposite",description="TODO",keywords={"vcf","disease","annotation","pedigree"})
public abstract class VCFComposite extends Launcher {
	private static final Logger LOG= Logger.build(VCFComposite.class).make();
	@Parameter(names={"-ped","--pedigree"},description="Pedigree file")
	private Pedigree pedigree=null;
	@Parameter(names={"-o","--out"},description="OUtput file")
	private File outputFile=null;
	
	private VCFFileReader vcfFileReader= null;
	
	protected List<Interval> intervalBuffer;
	protected void doWork(com.github.lindenb.jvarkit.util.vcf.VcfIterator r, VariantContextWriter w)
			throws IOException
			{
			VCFHeader header=r.getHeader();
			
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			//h2.addMetaDataLine(new VCFInfoHeaderLine(this.TAG,1,VCFHeaderLineType.String,"Values from bigwig file: "+BIGWIG));
			
			w.writeHeader(h2);

			while(r.hasNext())
				{
				
				r.next();
				}
			
		
			}
	@Override
	public int doWork(List<String> args) {
		
		try {
			this.vcfFileReader = new VCFFileReader(new File(oneAndOnlyOneFile(args)),true);
			
			this.vcfFileReader.close();
			this.vcfFileReader=null;
			return 0;
		} catch (Exception e) {
			return -1;
		} finally 
			{
			CloserUtil.close(vcfFileReader);
			}
		
		}
	
	}
