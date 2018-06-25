/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ReadNameSortMethod;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;

/** 

BEGIN_DOC

## Input

  * BAM : MUST be sorted using Picard SortSam (see https://github.com/samtools/hts-specs/issues/5 )
  * VCF : only SNP are considered. Genotypes are ignored (all ALT alleles are observed regardless of the sample/genotype )

##Example

```
$ java -jar picard.jar SortSam I=src/test/resources/S1.bam O=query.bam SO=queryname
$ java -jar dist/biostar322664.jar -V src/test/resources/S1.vcf.gz query.bam  


RF02_358_926_2:0:0_2:1:0_83	83	RF02	857	60	70M	=	358	-569	GACGTGAACTATATAATTAAAATGGACAGAAATCTGCCATCAACAGCTAGATATATAAGACCTAATTTAC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
RF02_362_917_2:0:0_2:1:0_6f	147	RF02	848	60	70M	=	362	-556	ATAAGGAATCACGTTAACTATATACTTAAAATGGACTGAAATCTGCCATCAACAGCTAGATATATAAGAC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
(...)
```
END_DOC

*/

@Program(name="biostar322664",
description="Extract PE Reads (with their mates) supporting variants in vcf file",
biostars=322664
)
public class Biostar322664 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar322664.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-V","--variant"},description="Variant VCF file. This tool **doesn't work** with INDEL/SV.",required=true)
	private File vcfFile = null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();

	
	
	@Override
	public int doWork(final List<String> args) {
		final IntervalTreeMap<VariantContext> variantMap=new IntervalTreeMap<>();
		VCFFileReader vcfFileReader=null;
		CloseableIterator<VariantContext> viter=null;
		SamReader samReader = null;
		SAMFileWriter samFileWriter = null;
		try {
			LOG.info("reading VCF "+this.vcfFile+" in memory...");
			vcfFileReader = new VCFFileReader(this.vcfFile, false);
			viter= vcfFileReader.iterator();
			viter.stream().
				filter(V->V.isVariant()).
				map(V->new VariantContextBuilder(V).attributes(Collections.emptyMap()).noGenotypes().make()).
				forEach(V->variantMap.put(new Interval(V.getContig(), V.getStart(), V.getEnd()),V));
			CloserUtil.close(viter);viter=null;
			vcfFileReader.close();vcfFileReader=null;
			LOG.info("Done Reading: "+this.vcfFile);

			
			samReader =super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader header = samReader.getFileHeader();
			if(header.getSortOrder()!=SAMFileHeader.SortOrder.queryname) {
				LOG.error("Expected SAM input to be sorted on "+SAMFileHeader.SortOrder.queryname+" but got "+header.getSortOrder());
				return -1;
				}
			samFileWriter = this.writingBamArgs.openSAMFileWriter(this.outputFile, header, true);
			CloseableIterator<SAMRecord> iter = samReader.iterator();
			EqualRangeIterator<SAMRecord> eq_range = new EqualRangeIterator<>(iter,ReadNameSortMethod.picard.get());
			while(eq_range.hasNext()) {
				boolean any_match=false;
				final List<SAMRecord> array= eq_range.next();
				
				for(final SAMRecord rec:array) {
					if(rec.getReadUnmappedFlag()) continue;
					final Collection<VariantContext> variants = variantMap.getOverlapping(rec);
					if(variants.isEmpty()) continue;
					final Cigar cigar=rec.getCigar();
					if(cigar==null || cigar.isEmpty()) continue;
					
					final byte bases[] = rec.getReadBases();
					if(bases==null || bases.length==0 || bases==SAMRecord.NULL_SEQUENCE) continue;
					
					for(final VariantContext ctx:variants) {
						if(ctx.getReference().length()!=1) continue;
						for(final Allele alt:ctx.getAlternateAlleles()) {
							if(alt.isNoCall() || alt.isSymbolic()) continue;
							if(alt.length()!=1) continue;
							int ref1 =rec.getUnclippedStart();
							int readpos = 0;
							
							for(final CigarElement ce:cigar) {
								if(ref1> ctx.getEnd()) break;
								final CigarOperator op =ce.getOperator();
								switch(op) {
								case P:break;
								case S: 
									{
									readpos+=ce.getLength();
									ref1+=ce.getLength();
									break;
									}
								case H:
									{
									ref1+=ce.getLength();
									break;
									}
								case I: readpos+=ce.getLength();break;
								case D:case N: ref1+=ce.getLength();break;
								case M:case EQ: case X: 
									{
									for(int x=0;x< ce.getLength();++x)
										{
										int rp2 = readpos+x;
										if(rp2>=0 && rp2<bases.length && ref1+x==ctx.getStart())
											{
											final char base = Character.toUpperCase((char)bases[rp2]);
											if(base==Character.toUpperCase(alt.getBases()[0]));
												{	
												System.err.println("Got "+ctx+" at "+rp2+" in "+rec.getContig()+" "+rec.getStart()+" "+rec.getEnd());
												any_match = true;
												break;
												}
											}
										}
									readpos+=ce.getLength();
									ref1+=ce.getLength();
									break;
									}
								default: throw new IllegalStateException("bad cigar operator "+op);
								}
							}
						if(any_match) break;
						}//end of cigar loop
					if(any_match) break;
					} // end of loop over variants
				}// end of for(Samrecord in array)
				if(any_match) {
					for(final SAMRecord rec:array) samFileWriter.addAlignment(rec);
					}
				
				}
			eq_range.close();
 			iter.close();
			samFileWriter.close();samFileWriter=null;
			samReader.close();samReader=null;
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
		} finally
			{
			CloserUtil.close(samReader);
			CloserUtil.close(samFileWriter);
			CloserUtil.close(vcfFileReader);
			}
	}
		
	
	public static void main(final String[] args) throws IOException
		{
		new Biostar322664().instanceMainWithExit(args);
		}
		

	}
