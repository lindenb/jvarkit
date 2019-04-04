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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfstripannot;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example

```
$ curl -sL "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
java -jar dist/vcfstripannot.jar -x "INFO/CSQ,INFO/EFF,INFO/AC,INFO/BaseQRankSum'  |\
grep -v "##"

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M10475	M10478	M10500	M128215
chr1	145273345	.	T	C	289.85	.	AF=0.38;DP=1000;DS;Dels=0.00;FS=3.974;HRun=1;HaplotypeScore=17.4275;MQ=29.25;MQ0=0;MQRankSum=-1.370;QD=0.39;ReadPosRankSum=-1.117	GT:AD:DP:GQ:PL	0/0:226,22:250:99:0,158,4259	0/1:224,24:250:5.77:6,0,5314	0/1:219,28:249:57.30:57,0,5027	0/1:215,34:250:99:269,0,3796
chr1	156011444	.	T	C	2523.46	.	AF=0.50;DP=204;Dels=0.00;FS=4.328;HRun=0;HaplotypeScore=4.3777;MQ=35.24;MQ0=0;MQRankSum=-0.101;QD=14.93;ReadPosRankSum=1.575	GT:AD:DP:GQ:PL	0/1:24,15:40:99:214,0,443	0/1:32,36:68:99:702,0,794	1/1:1,59:61:99:1656,132,0	0/0:34,1:35:69.10:0,69,717
chr5	64982321	.	T	C	61.12	.	AF=1.00;DP=4;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=0.0000;MQ=37.00;MQ0=0;QD=20.37	GT:AD:DP:GQ:PL	1/1:0,2:2:6:58,6,0	1/1:0,1:1:3.01:37,3,0	./.	./.
chr10	1142208	.	T	C	3404.30	.	AF=1.00;DP=122;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=2.6747;MQ=36.00;MQ0=0;QD=27.90	GT:AD:DP:GQ:PL	1/1:1,37:39:87.16:940,87,0	1/1:0,29:29:78.20:899,78,0	1/1:0,24:24:66.14:729,66,0	1/1:0,30:30:75.18:836,75,0
chr10	126678092	.	G	A	89.08	.	AF=0.13;DP=185;Dels=0.00;FS=3.490;HRun=0;HaplotypeScore=3.3843;MQ=25.32;MQ0=0;MQRankSum=6.568;QD=2.02;ReadPosRankSum=-5.871	GT:AD:DP:GQ:PL	0/0:64,3:67:99:0,165,1505	0/0:11,1:12:7.31:0,7,240	0/0:52,10:62:54.97:0,55,1263	0/1:35,9:44:99:125,0,693
chr10	135210791	.	T	C	65.41	.	AF=0.50;DP=11;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=0.2489;MQ=35.12;MQ0=0;MQRankSum=0.248;QD=16.35;ReadPosRankSum=-1.001	GT:AD:DP:GQ:PL	0/0:4,0:4:9:0,9,84	1/1:0,3:3:6.02:74,6,0	1/1:0,1:1:3.01:37,3,0	0/0:3,0:3:9.02:0,9,100
chr13	48873835	.	G	A	58.95	.	AF=1.00;DP=3;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=37.00;MQ0=0;QD=19.65	GT:AD:DP:GQ:PL	./.	./.	1/1:0,2:2:6.01:62,6,0	1/1:0,1:1:3.01:31,3,0
chr20	36779424	.	G	A	128.76	.	AF=0.13;DP=196;Dels=0.00;FS=1.447;HRun=0;HaplotypeScore=4.5749;MQ=36.22;MQ0=0;MQRankSum=-0.814;QD=3.90;ReadPosRankSum=-0.570	GT:AD:DP:GQ:PL	0/0:49,1:52:63.68:0,64,969	0/0:17,0:17:30.05:0,30,320	0/0:93,0:94:99:0,216,2384	0/1:24,9:33:99:165,0,505
chrX	17819377	.	T	C	7515.25	.	AF=1.00;DP=319;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=7.7850;MQ=36.33;MQ0=0;QD=23.56	GT:AD:DP:GQ:PL	1/1:0,125:126:99:2343,237,0	1/1:0,26:26:78.14:837,78,0	1/1:0,90:92:99:2640,244,0	1/1:0,74:75:99:1695,171,0
```

## History

* April 2017 : switched to BCFTOOLS syntax

END_DOC

 */
@Program(name="vcfstripannot",
	description="Removes one or more field from the INFO/FORMAT column of a VCF.",
	deprecatedMsg="Use bcftools annotate -x ",
	keywords={"vcf"}
	)
public class VCFStripAnnotations extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFStripAnnotations.class).make();	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	
	
	@XmlType(name="vcfstripannot")
	@XmlRootElement(name="vcfstripannot")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{
			@XmlElement(name="exclude")
			@Parameter(names={"-x","--exclude"},description="Use bcftools syntax INFO/x,INFO/y")
			private List<String> bcfToolStringSet=new ArrayList<>();
			
			private class CtxWriter extends DelegateVariantContextWriter
				{	
				private final Set<String> infoToKeep= new HashSet<>();
				private final Set<String> infoToRemove= new HashSet<>();
				private final Set<String> formatToKeep= new HashSet<>();
				private final Set<String> formatToRemove= new HashSet<>();
				private final Set<String> filterToKeep= new HashSet<>();
				private final Set<String> filterToRemove= new HashSet<>();
				private boolean removeId=false;
				
				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				
				@Override
				public void writeHeader(final VCFHeader header) {
					Set<VCFHeaderLine> vcfHeaderLines=  header.getMetaDataInInputOrder();
					for(String bcfToolString: CtxWriterFactory.this.bcfToolStringSet) {
						boolean inverse=false;
						if(bcfToolString.startsWith("^")) {
							inverse=true;
							bcfToolString=bcfToolString.substring(1);
							}
						for(final String bcfStr:bcfToolString.split("[,]"))
							{
							if(bcfStr.equals("ID"))
								{
								if(inverse)
									{
									LOG.warning("using inverse with ID");
									}
								this.removeId=true;
								continue;
								}
							else if(bcfStr.startsWith("INFO/"))
								{
								final String tag = bcfStr.substring(5);
								
								if(inverse)
									{
									this.infoToKeep.add(tag);
									}
								else
									{
									this.infoToRemove.add(tag);
									}
								}
							else if(bcfStr.startsWith("FILTER/"))
								{
								final String filter = bcfStr.substring(7);
								if(inverse)
									{
									this.filterToKeep.add(filter);
									}
								else
									{
									this.filterToRemove.add(filter);
									}
								}
							else if(bcfStr.startsWith("FORMAT/"))
								{
								final String format = bcfStr.substring(7);
								if(inverse)
									{
									this.formatToKeep.add(format);
									}
								else
									{
									this.formatToRemove.add(format);
									}
								continue;
								}
							else if(bcfStr.equals("FILTER"))
								{
								this.filterToRemove.add("*");
								}
							else if(bcfStr.equals("INFO"))
								{
								this.infoToRemove.add("*");
								}
							else
								{
								throw new JvarkitException.UserError("Cannot decode "+bcfStr+" in "+bcfToolString);
								}
							}
						}
					if( this.formatToKeep.contains(VCFConstants.GENOTYPE_KEY) ||
						this.formatToRemove.contains(VCFConstants.GENOTYPE_KEY))
						{
						throw new JvarkitException.UserError("Cannot remove/keep protected FORMAT:"+VCFConstants.GENOTYPE_KEY);
						}
					
					if(!this.filterToKeep.isEmpty() && !this.filterToRemove.isEmpty())
						{
						throw new JvarkitException.UserError("Cannot keep and remove FILTER at the same time");
						}
					if(!this.infoToKeep.isEmpty() && !this.infoToKeep.isEmpty())
						{
						throw new JvarkitException.UserError("Cannot keep and remove INFO at the same time");
						}

					if(!this.formatToKeep.isEmpty() && !this.formatToRemove.isEmpty())
						{
						throw new JvarkitException.UserError("Cannot keep and remove FORMAT at the same time");
						}
					if(!this.filterToKeep.isEmpty())
						{
						vcfHeaderLines = vcfHeaderLines.stream().
							filter(H->{
								if(!(H instanceof VCFFilterHeaderLine)) return true;
								final VCFFilterHeaderLine h = VCFFilterHeaderLine.class.cast(H);
								return filterToKeep.contains(h.getID());
								}).collect(Collectors.toSet());
						}
					if(!this.filterToRemove.isEmpty())
						{
						vcfHeaderLines = vcfHeaderLines.stream().
							filter(H->{
								if(!(H instanceof VCFFilterHeaderLine)) return true;
								if(filterToRemove.contains("*")) return false;
								final VCFFilterHeaderLine h = VCFFilterHeaderLine.class.cast(H);
								return !filterToRemove.contains(h.getID());
								}).collect(Collectors.toSet());
						}
					
					if(!this.infoToKeep.isEmpty())
						{
						vcfHeaderLines = vcfHeaderLines.stream().
							filter(H->{
								if(!(H instanceof VCFInfoHeaderLine)) return true;
								final VCFInfoHeaderLine h = VCFInfoHeaderLine.class.cast(H);
								return infoToKeep.contains(h.getID());
								}).collect(Collectors.toSet());
						}
					if(!this.infoToRemove.isEmpty())
						{
						vcfHeaderLines = vcfHeaderLines.stream().
							filter(H->{
								if(!(H instanceof VCFInfoHeaderLine)) return true;
								if(this.infoToRemove.contains("*")) return false;
								final VCFInfoHeaderLine h = VCFInfoHeaderLine.class.cast(H);
								return !this.infoToRemove.contains(h.getID());
								}).collect(Collectors.toSet());
						}

					if(!this.formatToKeep.isEmpty())
						{
						vcfHeaderLines = vcfHeaderLines.stream().
							filter(H->{
								if(!(H instanceof VCFInfoHeaderLine)) return true;
								final VCFFormatHeaderLine h = VCFFormatHeaderLine.class.cast(H);
								return this.formatToKeep.contains(h.getID());
								}).collect(Collectors.toSet());
						}
					if(!this.formatToRemove.isEmpty())
						{
						vcfHeaderLines = vcfHeaderLines.stream().
							filter(H->{
								if(!(H instanceof VCFFormatHeaderLine)) return true;
								final VCFFormatHeaderLine h = VCFFormatHeaderLine.class.cast(H);
								return !this.formatToRemove.contains(h.getID());
								}).collect(Collectors.toSet());
						}

					
					
					
					final VCFHeader h2= new VCFHeader(
							vcfHeaderLines,
							header.getSampleNamesInOrder()
							);
					//addMetaData(h2);
					
					super.writeHeader(h2);
					}
				
				@Override
				public void add(final VariantContext ctx) {
					final VariantContextBuilder b=new VariantContextBuilder(ctx);
					if(this.removeId) b.noID();
					
					/* INFO */
					if(!this.infoToKeep.isEmpty())
						{
						for(final String k2: ctx.getAttributes().keySet())
							{
							if(!this.infoToKeep.contains(k2)) b.rmAttribute(k2);
							}
						}
					if(!this.infoToRemove.isEmpty())
						{
						for(final String k2: ctx.getAttributes().keySet())
							{
							if(this.infoToRemove.contains(k2) || this.infoToRemove.contains("*")) {
								b.rmAttribute(k2);
								}
							}
						}
					/* formats */
					if(!this.formatToKeep.isEmpty())
						{
						final List<Genotype> genotypes=new ArrayList<Genotype>();
						for(final Genotype g:ctx.getGenotypes())
							{
							final GenotypeBuilder gb=new GenotypeBuilder(g);
							final Map<String, Object> map=new HashMap<String, Object>();
							for(final String key: g.getExtendedAttributes().keySet())
								{
								if(!this.formatToKeep.contains(key)) continue;
								map.put(key, g.getExtendedAttribute(key));
								}
							gb.attributes(map);
							genotypes.add(gb.make());
							}
						b.genotypes(genotypes);
						}
					if(!this.formatToRemove.isEmpty())
						{
						final List<Genotype> genotypes=new ArrayList<Genotype>();
						for(final Genotype g:ctx.getGenotypes())
							{
							final GenotypeBuilder gb=new GenotypeBuilder(g);
							if(this.formatToRemove.contains(VCFConstants.GENOTYPE_PL_KEY))
								{
								gb.noPL();
								}
							if(this.formatToRemove.contains("DP"))
								{
								gb.noDP();
								}
							if(this.formatToRemove.contains("AD"))
								{
								gb.noAD();
								}
							if(this.formatToRemove.contains("GQ"))
								{
								gb.noGQ();
								}
							if(this.formatToRemove.contains("*"))
								{
								gb.noAttributes();
								continue;
								}
							final Map<String, Object> map=new HashMap<String, Object>();
							for(final String key: g.getExtendedAttributes().keySet())
								{
								if(key.equals("DP") || key.equals("AD") || key.equals("GQ") || key.equals("PL")) continue;
								if(this.formatToRemove.contains(key)) continue;
								map.put(key, g.getExtendedAttribute(key));
								}
							gb.attributes(map);
							genotypes.add(gb.make());
							}
						b.genotypes(genotypes);
						}
					
					
					if(!this.filterToKeep.isEmpty())
						{
						final Set<String> filters = new HashSet<>( ctx.getFilters());
						filters.removeIf(S->!this.filterToKeep.contains(S));
						b.filters(filters);
						}
					if(!this.filterToRemove.isEmpty())
						{
						final Set<String> filters = new HashSet<>( ctx.getFilters());
						filters.removeAll(this.filterToRemove);
						if(this.filterToRemove.contains("*")) filters.clear();
						b.filters(filters);
						}

					super.add(b.make());
					}
				
				}
			
			@Override
			public VariantContextWriter open(final VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			
			}
	
	
	public VCFStripAnnotations()
		{
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter delegate
			) {	
		final VariantContextWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(iter.getHeader()).logger(LOG);
		out.writeHeader(iter.getHeader());
		while(iter.hasNext())
			{
			out.add(progress.watch(iter.next()));
			}
		out.close();
		progress.finish();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		
		}
	
	public static void main(String[] args) throws IOException
		{
		new VCFStripAnnotations().instanceMainWithExit(args);
		}
	}
