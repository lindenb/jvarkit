/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.tools.vcfcmp.EqualRangeVcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="vcfcalledwithanothermethod",description="After a VCF-merge, read a VCF, look back at some BAMS to tells if the missing genotypes were homozygotes-ref or not-called. If the number of reads is greater than min.depth, then the missing genotypes is said hom-ref.")
public class VcfCalledWithAnotherMethod extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfCalledWithAnotherMethod.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-f","--vcfs"},description="List of path to alternate VCF files")
	private File vcfList = null;
	//@ParametersDelegate
//	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	enum FlagGT{
		NO_FOUND_IN_SAMTOLS(0),
		SAME_IN_SAMTOLS(1),
		DISCORDANT_IN_SAMTOLS(-1)
		;
		
		int flag;
		FlagGT(int flag) {this.flag=flag;}
		
		static FlagGT get(int idx) {
			switch(idx) {
			case 0: return NO_FOUND_IN_SAMTOLS;
			case 1: return SAME_IN_SAMTOLS;
			case -1: return DISCORDANT_IN_SAMTOLS;
			default: throw new IllegalStateException();
			}
			}
		}
	
	enum FlagCtx{
		NO_FOUND_IN_SAMTOLS(0),
		SAME_IN_SAMTOLS(1),
		DISCORDANT_IN_SAMTOLS(-1)
		;
		
		int flag;
		FlagCtx(int flag) {this.flag=flag;}
		
		static  FlagCtx get(int idx) {
			switch(idx) {
			case 0: return NO_FOUND_IN_SAMTOLS;
			case 1: return SAME_IN_SAMTOLS;
			case -1: return DISCORDANT_IN_SAMTOLS;
			default: throw new IllegalStateException();
			}
			}

		}

	public VcfCalledWithAnotherMethod()
		{
		
		}
	@Override
	public int doWork(List<String> args) {
		if(this.vcfList==null) {
			LOG.error("Missing VCF");
			return -1;
		}
				
		SAMSequenceDictionary dictionary=null;
		VcfIterator in=null;
		VcfIterator samtoolsIterator=null;
		try {
			final Set<String> samtoolsFiles = 
					Files.readAllLines(this.vcfList.toPath()).stream().
					filter(L->!(L.startsWith("#") || L.trim().isEmpty())).
					collect(Collectors.toSet())
					;

			/** collect sequence dictionaries */
			for(final String samtoolsVcf: samtoolsFiles)
				{
				LOG.info("Reading header for "+samtoolsVcf);
				in  = VCFUtils.createVcfIterator(samtoolsVcf);
				final VCFHeader header= in.getHeader();
				final SAMSequenceDictionary dict = header.getSequenceDictionary();
				in.close();in=null;
				if(dict==null) 
					{
					LOG.error("vcf "+samtoolsVcf+" is missing dict");
					return -1;
					}
				if(dictionary == null)
					{
					dictionary=dict;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dictionary, dict))
					{
					LOG.error("not same sequence dictionary");
					return -1;
					}		
				}
			if(dictionary==null) {
				LOG.error("no dictionary, no input file.");
				return -1;
				}
			/** open primitive input */
			in =  super.openVcfIterator(oneFileOrNull(args));
			final StringBuilder description=new StringBuilder("Other Genotype generated by "+this.getProgramName()+":");
			for(FlagGT flag:FlagGT.values()) description.append(" "+flag.name()+"="+flag.flag);
			
			final VCFFormatHeaderLine formatHeaderLine = new VCFFormatHeaderLine(
					"IST"/* In SamTools */, 1,
					VCFHeaderLineType.Integer,
					description.toString()
					);
			
			description.setLength(0);
			description.append("Other VariantContext generated by "+this.getProgramName()+":");
			for(FlagCtx flag:FlagCtx.values()) description.append(" "+flag.name()+"="+flag.flag);
			
			final VCFInfoHeaderLine infoHeaderLine = new VCFInfoHeaderLine(
					formatHeaderLine.getID(),1,
					VCFHeaderLineType.Integer,
					description.toString()
					);
			description.setLength(0);
			final File tmpDir= new File(System.getProperty("java.io.tmpdir"));
			final File tmpFile1 = File.createTempFile("fixvcf", ".vcf",tmpDir);
			final File tmpFile2 = File.createTempFile("fixvcf", ".vcf",tmpDir);
			final VCFHeader header=in.getHeader();
			VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(formatHeaderLine);
			h2.addMetaDataLine(infoHeaderLine);
			super.addMetaData(h2);
			
			int loopIndex=0;
			for(String samtoolsVcf: samtoolsFiles)
				{				
				
				LOG.info("VCF: "+samtoolsVcf);
				samtoolsIterator = VCFUtils.createVcfIterator(samtoolsVcf);
				final EqualRangeVcfIterator equalRange = new EqualRangeVcfIterator(samtoolsIterator,
						VCFUtils.createTidPosRefComparator(dictionary)
						);

				final VariantContextWriter w = VCFUtils.createVariantContextWriter(loopIndex%2==0?tmpFile1:tmpFile2);
				
				w.writeHeader(h2);
				while(in.hasNext())
					{
					final VariantContext ctx= in.next();
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);

					FlagCtx ctxflag = FlagCtx.get((Integer)ctx.getAttributeAsInt(infoHeaderLine.getID(),FlagCtx.NO_FOUND_IN_SAMTOLS.flag));

					final List<VariantContext> samtoolsCtxList = equalRange.next(ctx);
					
					/* no sample in main vcf */
					if(!h2.hasGenotypingData())
						{
						/* just check same ALT found */
						if(!samtoolsCtxList.isEmpty() &&  ctxflag!=FlagCtx.DISCORDANT_IN_SAMTOLS) {
							final Set<Allele> alt1s = new TreeSet<>(ctx.getAlternateAlleles());
							final Set<Allele> alt2s = new TreeSet<>();
							for(final VariantContext samtoolsCtx:samtoolsCtxList)
								{
								alt2s.addAll(samtoolsCtx.getAlternateAlleles());
								}
							if(!alt1s.equals(alt2s)) {
								ctxflag = FlagCtx.DISCORDANT_IN_SAMTOLS;
								}
							else
								{
								ctxflag = FlagCtx.SAME_IN_SAMTOLS;
								}
							}
						vcb.attribute(infoHeaderLine.getID(), ctxflag.flag);
						w.add(vcb.make());
						continue;
						}
					
					
					final List<Genotype> genotypes = new ArrayList<>(header.getNGenotypeSamples());
					
					for(final String sample:header.getSampleNamesInOrder()) {
						final Genotype gatkGenotype = ctx.getGenotype(sample);
						Object ogt = gatkGenotype.getExtendedAttribute(formatHeaderLine.getID(),FlagGT.NO_FOUND_IN_SAMTOLS.flag);
						if(ogt instanceof String)
							{
							ogt =Integer.parseInt(ogt.toString());
							}
						
						
						FlagGT gtflag =  FlagGT.get((Integer)ogt);
						
						final GenotypeBuilder gb = new GenotypeBuilder(gatkGenotype);
						for(final VariantContext samtoolsCtx:samtoolsCtxList) {
							final Genotype stGenotype = samtoolsCtx.getGenotype(sample);
							
							if(stGenotype!=null) //this samtools vcf contains the sample
								{
								boolean sameGt = gatkGenotype.sameGenotype(stGenotype,true/*ignore phase*/);
								if(sameGt && gtflag!=FlagGT.DISCORDANT_IN_SAMTOLS) {
									gtflag= FlagGT.SAME_IN_SAMTOLS;
									}
								else if(!sameGt)
									{
									gtflag= FlagGT.DISCORDANT_IN_SAMTOLS;
									}
								/* update context flag */
								if(gtflag==FlagGT.DISCORDANT_IN_SAMTOLS) 
									{
									ctxflag = FlagCtx.DISCORDANT_IN_SAMTOLS;
									}
								else if(ctxflag!=FlagCtx.DISCORDANT_IN_SAMTOLS &&
										gtflag == FlagGT.SAME_IN_SAMTOLS)
									{
									ctxflag = FlagCtx.SAME_IN_SAMTOLS;
									}
								
								}
							}
						
						gb.attribute(formatHeaderLine.getID(), gtflag.flag);
						genotypes.add(gb.make());
						}
					
					
				
					vcb.attribute(infoHeaderLine.getID(), ctxflag.flag);
					vcb.genotypes(genotypes);
					w.add(vcb.make());
					}
				w.close();
				in.close();
				equalRange.close();
				CloserUtil.close(samtoolsIterator);
				
				LOG.info("done vcf "+samtoolsVcf);
				
				//reopen in
				in = VCFUtils.createVcfIteratorFromFile(loopIndex%2==0?tmpFile1:tmpFile2);
				++loopIndex;
				h2= in.getHeader();
				}
			
			final VariantContextWriter w = super.openVariantContextWriter(null,outputFile);
			w.writeHeader(h2);
			while(in.hasNext())
				{
				w.add(in.next());
				}
			in.close();
			w.close();
			tmpFile1.delete();
			tmpFile2.delete();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	public static void main(String[] args) {
		new VcfCalledWithAnotherMethod().instanceMainWithExit(args);

	}

}
