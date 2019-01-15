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
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
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
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;


/**
BEGIN_DOC
## Example

```bash
$  java -jar dist/vcfcomm.jar < in.vcf > out.vcf
``

END_DOC
 */
@Program(name="vcfcomm",description="Equivalent of linux comm for VCF")
public class VCFComm extends AbstractVCFCompareBase {
	private final Logger LOG=Logger.build(VCFComm.class).make();
	private final Pattern tab = Pattern.compile("[\t]");
	public VCFComm() 
		{
		}
	@Parameter(names="-a",description="ignore variations present in ALL files")
	private boolean ignore_everywhere=false;
	@Parameter(names="-A",description="only print variations present in ALL files")
	private boolean only_everywhere=false;
	@Parameter(names={"-norm","--normalize"},description="normalize chromosomes names (remove chr prefix, chrM -> MT)")
	private boolean normalize_chr=false;

	
	/* we can remove INFO from the line */
	@Override
	protected String simplify(final String line, int input_idx) {
		final String tokens[]=tab.split(line);
		if(normalize_chr){
			String chr=tokens[0];
			if(chr.startsWith("chr")) chr=chr.substring(3);
			if(chr.equals("M")) chr="MT";
			tokens[0]=chr;
			}
		if(tokens.length>7)
			{
			tokens[7]=VCFConstants.MISSING_VALUE_v4;
			return String.join(VCFConstants.FIELD_SEPARATOR, tokens);
			}
		else
			{
			return line;
			}	
		}
	
	@Override
	public int doWork(final List<String> args) {
		CloseableIterator<LineAndFile> iter = null;
		SortingCollection<LineAndFile> variants=null;
		VariantContextWriter w=null;
		try
			{
			if(args.isEmpty())
				{
				LOG.error("Illegal number of arguments");
				return -1;
				}
			
			Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
			
			
			
			variants=SortingCollection.newInstance(
					LineAndFile.class, 
					new LineAndFileCodec(),
					new LineAndFileComparator(),
					super.sortingCollectionArgs.getMaxRecordsInRam(),
					super.sortingCollectionArgs.getTmpPaths()
					);
			variants.setDestructiveIteration(true);
			
			/** new sample names in the output vcf: one  sample per file */
			final Map<Integer,String> fileid2sampleName=new TreeMap<>();
			/** samples names as they appear in the original VCF headers*/
			final Counter<String> countInputSamples=new Counter<String>();
		
			/** dicts */
			final List<SAMSequenceDictionary> all_dictionaries=new ArrayList<>();
			
			for(final String vcffilename:IOUtils.unrollFiles(args))
				{
				LOG.info("Reading from "+vcffilename);
				final Input input=super.put(variants, vcffilename);
				
				String sampleName=vcffilename;
				if(sampleName.endsWith(".vcf.gz"))
					{
					sampleName = sampleName.substring(0, sampleName.length()-7);
					}
				else if(sampleName.endsWith(".vcf.gz"))
					{
					sampleName = sampleName.substring(0, sampleName.length()-4);
					}
				int slash=sampleName.lastIndexOf(File.separatorChar);
				if(slash!=-1) sampleName=sampleName.substring(slash+1);
				int suffix=1;
				// loop until we find a uniq name
				for(;;)
					{
					final String key=sampleName+(suffix==1?"":"_"+suffix);
					if(fileid2sampleName.values().contains(key))
						{
						suffix++;
						continue;
						}
					fileid2sampleName.put(input.file_id, key);
					metaData.add(new VCFHeaderLine(key,vcffilename));
					break;
					}
				
				for(final String sname:input.codecAndHeader.header.getSampleNamesInOrder())
					{
					countInputSamples.incr(sname);
					}
				all_dictionaries.add(input.codecAndHeader.header.getSequenceDictionary());
				}
			
			
			
			
			variants.doneAdding();
			
			/** unique sample name, if any present in all VCF*/
			Optional<String> unqueSampleName=Optional.empty();
			if(countInputSamples.getCountCategories()==1 &&
				countInputSamples.count(countInputSamples.keySet().iterator().next())==fileid2sampleName.size())
				{
				unqueSampleName=Optional.of(countInputSamples.keySet().iterator().next());
				LOG.info("Unique sample name is "+unqueSampleName.get());
				}
			
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.GENOTYPE_QUALITY_KEY,
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.GENOTYPE_FILTER_KEY)
					;
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.ALLELE_COUNT_KEY,
					VCFConstants.ALLELE_NUMBER_KEY
					);
			
			metaData.add(new VCFHeaderLine(getClass().getSimpleName(),"version:"+getVersion()+" command:"+getProgramCommandLine()));

			final VCFFilterHeaderLine variantNotCalledInAllVcf=new VCFFilterHeaderLine("NotCalledEveryWhere", "Variant was NOT called in all input VCF");
			metaData.add(variantNotCalledInAllVcf);
			final VCFFilterHeaderLine variantWasFiltered = new VCFFilterHeaderLine("VariantWasFiltered", "At least one variant was filtered");
			metaData.add(variantWasFiltered);
			final VCFFormatHeaderLine variantQUALFormat = new VCFFormatHeaderLine("VCQUAL",1,VCFHeaderLineType.Float,"Variant Quality");
			metaData.add(variantQUALFormat);
			metaData.add(new VCFFormatHeaderLine(VCFConstants.ALLELE_NUMBER_KEY,1,VCFHeaderLineType.Integer,"Number of allle in the src vcf"));
			metaData.add(new VCFFormatHeaderLine(VCFConstants.ALLELE_COUNT_KEY,1,VCFHeaderLineType.Integer,"Number of ALT alllele"));
			
			
			final VCFInfoHeaderLine foundInCountVcfInfo = new VCFInfoHeaderLine("NVCF",1,VCFHeaderLineType.Integer,"Number of VCF this variant was found");
			metaData.add(foundInCountVcfInfo);
			final VCFInfoHeaderLine variantTypesInfo = new VCFInfoHeaderLine("VTYPES",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"Distinct Variants type");
			metaData.add(variantTypesInfo);
			final VCFFilterHeaderLine multipleTypeFilters = new VCFFilterHeaderLine("DiscordantTypes", "Discordant types at this position");
			metaData.add(multipleTypeFilters);
			final VCFFormatHeaderLine variantTypeFormat = new VCFFormatHeaderLine("VTYPE",1,VCFHeaderLineType.String,"Variant Type");
			metaData.add(variantTypeFormat);
			
			
			final VCFFilterHeaderLine uniqueVariantDiscordantGTFilter;

			
			if(unqueSampleName.isPresent())
				{
				metaData.add(new VCFHeaderLine("UniqSample",unqueSampleName.get()));
				uniqueVariantDiscordantGTFilter = new VCFFilterHeaderLine("DiscordantGenotypeForUniqSample", "Genotype Dicordant for for sample "+unqueSampleName.get());
				metaData.add(uniqueVariantDiscordantGTFilter);
				}
			else
				{
				uniqueVariantDiscordantGTFilter = null;
				}
			
			
			final VCFHeader header=new VCFHeader(
					metaData,
					new ArrayList<>(fileid2sampleName.values())
					);

			
			if(!normalize_chr && !all_dictionaries.contains(null))//all have a dict
				{
				SAMSequenceDictionary thedict=null;
				for(int x=0;x< all_dictionaries.size();++x)
					{
					SAMSequenceDictionary d=all_dictionaries.get(x);
					if(thedict==null )
						{
						thedict=d;
						}
					else if(!SequenceUtil.areSequenceDictionariesEqual(d, thedict))
						{
						thedict=null;
						break;
						}
					}
				if(thedict!=null) header.setSequenceDictionary(thedict);
				}
			
			
			 w= super.openVariantContextWriter(super.outputFile);
			w.writeHeader(header);
			final List<LineAndFile> row=new ArrayList<LineAndFile>(super.inputs.size());
			
			final Comparator<LineAndFile> posCompare = (A,B)->A.getContigPosRef().compareTo(B.getContigPosRef());

			
			
			iter=variants.iterator();
			
			for(;;)
				{
				LineAndFile rec=null;
				if(iter.hasNext())
					{
					rec=iter.next();
					}
			
				if(rec==null || (!row.isEmpty() && posCompare.compare(row.get(0),rec)!=0))
					{
					if(!row.isEmpty())
						{
						final VariantContext first=row.get(0).getContext();
						/* in which file id we find this variant */
						Set<Integer> fileids_for_variant= row.stream().map(LAF->LAF.fileIdx).collect(Collectors.toSet());
						
						// see with HAS multiple chrom/pos/ref but different alt
						if(row.size()!=fileids_for_variant.size())
							{
							for(;;)
								{
								boolean ok=true;
								for(int x=0;ok && x+1<row.size();++x)
									{
									final VariantContext ctxx=row.get(x).getContext();
									final List<Allele> altsx=ctxx.getAlternateAlleles();
									for(int y=x+1;ok && y< row.size();++y)
										{
										if(row.get(x).fileIdx!=row.get(y).fileIdx) continue;
										final VariantContext ctxy=row.get(y).getContext();
										final List<Allele> altsy=ctxy.getAlternateAlleles();
										if(altsx.equals(altsy)) continue;
										if(!ctxx.isVariant() && ctxy.isVariant())
											{
											row.remove(x);
											}
										else if(ctxx.isVariant() && !ctxy.isVariant())
											{
											row.remove(y);
											}
										else if(!ctxx.isSNP() && ctxy.isSNP())
											{
											row.remove(x);
											}
										else if(ctxx.isSNP() && !ctxy.isSNP())
											{
											row.remove(y);
											}
										else if(altsx.size()> altsy.size()) {
											row.remove(x);
											}
										else if(altsx.size()< altsy.size()) {
											row.remove(y);
											}
										else
											{
											row.remove(y);
											}
										ok=false;
										break;
										}
									}
								if(ok) break;
								}
							fileids_for_variant= row.stream().
									map(LAF->LAF.fileIdx).
									collect(Collectors.toSet());
							}	
						if(row.size()!=fileids_for_variant.size())
							{
							LOG.error("There are some duplicated variants at the position "+new ContigPosRef(first)+" in the same vcf file");
							for(final LineAndFile laf:row)
								{
								LOG.error("File ["+laf.fileIdx+"]"+fileid2sampleName.get(laf.fileIdx));
								LOG.error("\t"+laf.getContigPosRef());
								}
							row.clear();
							}
						else
							{
							final Set<Allele> alleles= row.stream().
									flatMap(R->R.getContext().getAlleles().
									stream()).collect(Collectors.toSet());
							
							final VariantContextBuilder vcb=new VariantContextBuilder(
									getClass().getName(),
									first.getContig(),
									first.getStart(),
									first.getEnd(),
									alleles
									);
							final Set<String> filters = new HashSet<>();
							final Set<VariantContext.Type> variantContextTypes=new HashSet<>();
							final List<Genotype> genotypes=new ArrayList<Genotype>();
							for(final LineAndFile laf:row)
								{
								if(laf.getContext().isFiltered()) filters.add(variantWasFiltered.getID());
								variantContextTypes.add(laf.getContext().getType());
								final GenotypeBuilder gbuilder=new GenotypeBuilder();
								gbuilder.name(fileid2sampleName.get(laf.fileIdx));
								if(unqueSampleName.isPresent())
									{
									final Genotype g0=laf.getContext().getGenotype(unqueSampleName.get());
									if(g0==null){
										iter.close();
										w.close();
										throw new IllegalStateException("Cannot find genotype for "+unqueSampleName.get());
									}
									if(g0.hasDP()) gbuilder.DP(g0.getDP());
									if(g0.hasGQ()) gbuilder.GQ(g0.getGQ());
									gbuilder.alleles(g0.getAlleles());
									}
								else
									{
									gbuilder.alleles(Arrays.asList(first.getReference(),first.getReference()));
									if(laf.getContext().hasAttribute(VCFConstants.DEPTH_KEY))
										{
										gbuilder.DP(laf.getContext().getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
										}
									}
								if(laf.getContext().isFiltered() )
									{
									gbuilder.filter("VCFFILTERED");
									}
								if(laf.getContext().hasLog10PError())
									{
									gbuilder.attribute(variantQUALFormat.getID(), laf.getContext().getPhredScaledQual());
									}
								
								gbuilder.attribute(VCFConstants.ALLELE_NUMBER_KEY, 
										laf.getContext().getGenotypes().stream().flatMap(G->G.getAlleles().stream()).filter(A->!A.isNoCall()).count()
										);
								gbuilder.attribute(VCFConstants.ALLELE_COUNT_KEY, 
										laf.getContext().getGenotypes().stream().flatMap(G->G.getAlleles().stream()).filter(A->!(A.isReference() || A.isNoCall())).count()
										);
								gbuilder.attribute(variantTypeFormat.getID(),laf.getContext().getType().name());
								
								genotypes.add(gbuilder.make());
								}
							final String id=String.join(";",row.stream().map(LAF->LAF.getContext()).
									filter(V->V.hasID()).map(V->V.getID()).
									collect(Collectors.toSet()))
									;
							if(!id.isEmpty()) vcb.id(id);
							
							
							vcb.genotypes(genotypes);
							
							if(unqueSampleName.isPresent())
								{
								boolean all_same=true;
								for(int x=0;all_same && x+1< genotypes.size();++x)
									{
									if(!genotypes.get(x).isCalled()) continue;
									for(int y=x+1;all_same && y< genotypes.size();++y)
										{
										if(!genotypes.get(y).isCalled()) continue;
										if(!genotypes.get(x).sameGenotype(genotypes.get(y),true))
											{
											all_same=false;
											break;
											}
										}
									}
								if(!all_same) filters.add(uniqueVariantDiscordantGTFilter.getID());
								}
							
							//Add AN
							vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,
									genotypes.stream().
										filter(G->G.isCalled()).
										mapToInt(G->G.getAlleles().size()).sum()
									);
							if(!variantContextTypes.isEmpty()) {
								vcb.attribute(variantTypesInfo.getID(),
									new ArrayList<>(variantContextTypes.stream().map(T->T.name()).collect(Collectors.toSet()))
									);
								if(variantContextTypes.size()>1)
									{
									filters.add(multipleTypeFilters.getID());
									}
								}
							vcb.attribute(foundInCountVcfInfo.getID(), fileids_for_variant.size());
							
							boolean print=true;
							if(row.size()==super.inputs.size() && ignore_everywhere)
								{
								print=false;
								}
							if(fileids_for_variant.size()!=fileid2sampleName.size())
								{
								filters.add(variantNotCalledInAllVcf.getID());
								if(only_everywhere)
									{
									print=false;
									}
								}
							
							vcb.filters(filters);
							
							if(print)
								{
								w.add(vcb.make());
								}
							}
						row.clear();
						}
					if(rec==null) break;
					}
				
				row.add(rec);
				}
			iter.close();
			iter=null;
			
			w.close();
			w=null;

			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(w);
			try
				{
				if(variants!=null) variants.cleanup();
				}
			catch(Exception err)
				{
				}
			}
		}
	public static void main(String[] args)
		{
		new VCFComm().instanceMainWithExit(args);
		}
	}
