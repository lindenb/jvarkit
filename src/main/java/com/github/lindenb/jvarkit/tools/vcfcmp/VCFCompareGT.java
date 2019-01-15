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


*/
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BiFunction;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;

/**
BEGIN_DOC

## Input

input is a set of VCF files or a file with '.list' suffix with the path (one path per line).

Genotypes are supposed diploids.

## Example

```bash

$ java -jar dist/vcfcomparegt.jar -m  Sample.samtools.vcf.gz Sample.gatk.vcf.gz

##fileformat=VCFv4.1
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GCH,Number=1,Type=Integer,Description="Changed Genotype">
##FORMAT=<ID=GNW,Number=1,Type=Integer,Description="Genotype Created/Deleted">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Qual">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=GDF,Number=.,Type=String,Description="Samples with Genotype Difference">
##VCFCompareGT_1=File: Sample.samtools.vcf.gz
##VCFCompareGT_2=File: Sample.gatk.vcf.gz
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_2	Sample_1
X	1860854	rs5781	A	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:2:0:1:6	./.
X	1866893	rs2824	G	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:2:0:1:6	./.
X	1878904	.	G	C	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	0/1:20:0:1:71	./.
X	1895117	.	A	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/0:2:0:1:27
X	1895755	.	C	AG	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:4:0:1:17
X	1900009	rs6181	A	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:13:0:1:30	./.
X	1905130	.	AG	A	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:3:0:1:16
X	1905160	.	A	T	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:3
X	1905165	.	C	G	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:4
X	1913889	.	C	A	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:3
X	1948846	rs6	T	TG	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	1/1:239:0:1:99	./.
X	1955199	.	C	T	.	.	GDF=Sample	GT:DP:GCH:GNW:GQ	./.	1/1:1:0:1:4
(...)
```

END_DOC
 */
@Program(
	name="vcfcomparegt",
	description=" compare two or more genotype-callers for the same individuals. Produce a VCF with FORMAT fields indicating if a genotype is new or modified.",
	keywords={"vcf","compare"}
	)
public class VCFCompareGT extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFCompareGT.class).make();
	@Parameter(names="-m",description="only print modified samples")
	private boolean only_print_modified=false;

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-label","--labels"},description="A comma separated list of label that will be used as the title of the vcfs. Must be provided in the same order. If blank, some numeric indexes will be used")
	private String vcfLabelsStr ="";
	@Parameter(names={"-vf","--variant-filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = JexlVariantPredicate.create("");
	@Parameter(names={"-nc","--nocall2homref"},description="convert no call to hom-ref")
	private boolean convertNoCallToHomRef=false;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	private class Variant
		{
		String chrom="";
		String id=VCFConstants.EMPTY_ID_FIELD;
		String ref=VCFConstants.EMPTY_ALLELE;
		int start=-1;
		int end=-1;
		int file_index1 = 0;
		String sampleName="";
		String a1=VCFConstants.EMPTY_ALLELE;
		String a2=VCFConstants.EMPTY_ALLELE;
		int dp=-1;
		int gq=-1;
		boolean isHomRef() {
			return a1.equals(ref) && a1.equals(a2);
			}
		}
	
	private class PosComparator implements Comparator<Variant>
		{
		@Override
		public int compare(final Variant v1,final  Variant v2) {
			int i=v1.chrom.compareTo(v2.chrom);
			if(i!=0) return i;
			i=v1.start - v2.start;
			if(i!=0) return i;
			i=v1.end - v2.end;
			if(i!=0) return i;
			i=v1.ref.compareToIgnoreCase(v2.ref);
			if(i!=0) return i;
			return 0;
			}
		}
	
	private class VariantComparator extends PosComparator
		{
		@Override
		public int compare(final Variant v1,final  Variant v2) {
			int i=super.compare(v1, v2);
			if(i!=0) return i;
			i=v1.sampleName.compareTo(v2.sampleName);
			if(i!=0) return i;
			i= Integer.compare(v1.file_index1 ,v2.file_index1);
			if(i!=0) return i;
			return 0;
			}
		}
	
	private class VariantCodec extends AbstractDataCodec<Variant>
		{
		@Override
		public Variant decode(final DataInputStream dis) throws IOException
			{
			final Variant v=new Variant();
			try {
				v.chrom=dis.readUTF();
			} catch (final Exception e) {
				return null;
				}
			v.start=dis.readInt();
			v.end=dis.readInt();
			v.id=dis.readUTF();
			v.ref=dis.readUTF();
			v.sampleName=dis.readUTF();
			v.a1=readString(dis);
			v.a2=readString(dis);
			v.file_index1 = dis.readInt();
			v.dp=dis.readInt();
			v.gq=dis.readInt();
			return v;
			}
		@Override
		public void encode(final DataOutputStream dos,final Variant v)
				throws IOException
			{
			dos.writeUTF(v.chrom);
			dos.writeInt(v.start);
			dos.writeInt(v.end);
			dos.writeUTF(v.id);
			dos.writeUTF(v.ref);
			dos.writeUTF(v.sampleName);
			writeString(dos,v.a1);
			writeString(dos,v.a2);
			dos.writeInt(v.file_index1);
			dos.writeInt(v.dp);
			dos.writeInt(v.gq);
			}
		@Override
		public AbstractDataCodec<Variant> clone() {
			return new VariantCodec();
			}
		}
	
	
	public VCFCompareGT()
		{
		}
	

	@Override
	public int doWork(final List<String> arguments) {
		
		if(arguments.isEmpty())
			{
			LOG.error("VCFs missing.");
			return -1;
			}
		final List<String> vcfLabelList;
		if(!StringUtil.isBlank(this.vcfLabelsStr)) {
			vcfLabelList = Arrays.stream(this.vcfLabelsStr.split("[,]")).filter(S->!StringUtil.isBlank(S)).collect(Collectors.toList());
			if(new HashSet<>(vcfLabelList).size()!=arguments.size()) {
				LOG.error("bad number of labels in : "+this.vcfLabelsStr);
				return -1;
				}
			if(vcfLabelList.stream().anyMatch(S->!S.matches("[0-9A-Za-z]+")))
				{
				LOG.error("bad label in : "+this.vcfLabelsStr);
				return -1;
				}
			}
		else
			{
			vcfLabelList = new ArrayList<>();
			for(int i=0;i< arguments.size();++i) {
				vcfLabelList.add(String.valueOf(i+1));
				}
			}
		
		VariantComparator varcmp=new VariantComparator();
		SortingCollection<Variant> variants = null;
		final Set<String> sampleNames=new LinkedHashSet<>();
		try
			{
			
			variants=SortingCollection.newInstance(
					Variant.class,
					new VariantCodec(),
					varcmp,
					writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths()
					);
			variants.setDestructiveIteration(true);
			
			final Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
			metaData.add(new VCFHeaderLine(getClass().getSimpleName(),"version:"+getVersion()+" command:"+getProgramCommandLine()));
			final BiFunction<String, Integer, String> createName=(SN,IDX) -> 
				SN+"_"+vcfLabelList.get(IDX)
			;
			
			for(int i=0;i< arguments.size();++i)
				{
				final File vcfFile= new File(arguments.get(i));
				
				LOG.info("Opening "+vcfFile);
				final VCFFileReader vcfFileReader = new VCFFileReader(vcfFile,false);
				final CloseableIterator<VariantContext> iter = vcfFileReader.iterator();
				final VCFHeader header = vcfFileReader.getFileHeader();
				sampleNames.addAll(header.getSampleNamesInOrder());
				
				metaData.add(new VCFHeaderLine(
						getClass().getSimpleName()+"_"+vcfLabelList.get(i),
						"File: "+vcfFile.getPath())
						);

				
				long nLines=0;
				while(iter.hasNext())
					{
					final VariantContext var = iter.next();
					
					if(nLines++%10000==0)
						{
						LOG.info(vcfFile+" "+nLines);
						}
					if(!this.variantFilter.test(var)) continue;
					if(!var.isVariant()) continue;
					if(!var.hasGenotypes()) continue;
					for(final Genotype genotype:var.getGenotypes())
						{
						final Variant rec=new Variant();
						
						rec.file_index1= i+1;
						rec.sampleName=genotype.getSampleName();
						rec.chrom = var.getContig();
						rec.start = var.getStart();
						rec.end = var.getEnd();
						
						rec.ref = var.getReference().getDisplayString().toUpperCase();
						if(var.hasID())
							{
							rec.id=var.getID();
							}
						if(genotype.hasDP())
							{
							rec.dp=genotype.getDP();
							}
						if(genotype.hasGQ())
							{
							rec.gq=genotype.getGQ();
							}	
						final List<Allele> alleles=genotype.getAlleles();

						if(genotype.isNoCall() || !genotype.isAvailable() || alleles==null) {
							if(!this.convertNoCallToHomRef) continue;
							rec.a1 = var.getReference().getDisplayString().toUpperCase();
							rec.a2 = rec.a1;
							}
						else if(alleles.size()==1)
							{
							rec.a1=alleles.get(0).getDisplayString().toUpperCase();
							rec.a2=rec.a1;
							}
						else if(alleles.size()==2)
							{
							rec.a1=alleles.get(0).getDisplayString().toUpperCase();
							rec.a2=alleles.get(1).getDisplayString().toUpperCase();
							if(rec.a1.compareTo(rec.a2)>0)
								{
								final String tmp=rec.a2;
								rec.a2=rec.a1;
								rec.a1=tmp;
								}
							}
						else
							{
							continue;
							}
						variants.add(rec);
						}
					}
				iter.close();
				vcfFileReader.close();
				}
			variants.doneAdding();
		
			LOG.info("Done Adding");
			
			
			final String GenpotypeChangedKey="GCH";
			final String GenpotypeCreated="GNW";
			final String GenpotypeDiff="GDF";
			
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY));
			
			metaData.add(new VCFFormatHeaderLine(GenpotypeChangedKey,1,VCFHeaderLineType.Integer, "Changed Genotype"));
			metaData.add(new VCFFormatHeaderLine(GenpotypeCreated,1,VCFHeaderLineType.Integer, "Genotype Created/Deleted"));
			metaData.add(new VCFInfoHeaderLine(GenpotypeDiff,VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String, "Samples with Genotype Difference"));
			metaData.add(new VCFFilterHeaderLine("DISCORDANCE","something has changed."));

            final Set<String> newSampleNames=new TreeSet<>();
            for(int i=0;i< vcfLabelList.size();++i)
                   {
                   for(final String sample:sampleNames)
                           {
                           newSampleNames.add(createName.apply(sample,i));
                           }
                   }
			
			final VCFHeader header=new VCFHeader(
					metaData,
					new ArrayList<>(newSampleNames)
					);
			
			final VariantContextWriter w= super.openVariantContextWriter(outputFile);
			w.writeHeader(header);
			
			final PosComparator posCompare=new PosComparator();
			final EqualRangeIterator<Variant> iter=new EqualRangeIterator<>(variants.iterator(),posCompare);
			while(iter.hasNext())
				{
				final List<Variant> row= iter.next();
				/** this sample is not always the same */
				final Set<String> samplesModified=new TreeSet<>();
				/** the number of sample is different from vcflist.size() */
				final Set<String> samplesCreates=new TreeSet<>();
				
				final Map<String,List<Variant>> sample2variants = row.stream().collect(Collectors.groupingBy(T->T.sampleName));
				for(final String sn:sample2variants.keySet())
					{
					boolean all_hom_ref = true;
					final List<Variant> sampleVariants = sample2variants.get(sn);

					for(int x=0;x /*+1 non, besoin de tester hom_ref */ < sampleVariants.size();++x)
						{
						final Variant var1=sampleVariants.get(x);
						if(!var1.isHomRef()) all_hom_ref = false;
						for(int y=x+1;y< sampleVariants.size();++y)
							{
							final Variant var2=sampleVariants.get(y);
							if(var1.a1.equals(var2.a1) && var1.a2.equals(var2.a2) ) continue;
							samplesModified.add(var1.sampleName);
							}
						}
					
					if(sampleVariants.size() != arguments.size())
						{
						if(!convertNoCallToHomRef || (this.convertNoCallToHomRef && !all_hom_ref))
							{
							samplesCreates.add(sn);
							}
						}					
					}
				
				
				final Variant first = row.get(0);
				final Set<Allele> alleles = new HashSet<>();
				alleles.add(Allele.create(first.ref, true));
				for(final Variant var:row)
					{
					alleles.add(Allele.create(var.a1, var.a1.equalsIgnoreCase(var.ref)));
					alleles.add(Allele.create(var.a2, var.a2.equalsIgnoreCase(var.ref)));
					}
				
				
				final VariantContextBuilder b = new VariantContextBuilder(
						getClass().getName(),
						first.chrom,
						first.start,
						first.end,
						alleles
						);
				
				//build genotypes
				final List<Genotype> genotypes=new ArrayList<Genotype>();
				for(final Variant var:row)
					{
					//alleles for this genotype
					final List<Allele> galleles=new ArrayList<Allele>();
					galleles.add(Allele.create(var.a1, var.a1.equalsIgnoreCase(var.ref)));
					galleles.add(Allele.create(var.a2, var.a2.equalsIgnoreCase(var.ref)));
					
					final GenotypeBuilder gb=new GenotypeBuilder();
					gb.DP(var.dp);
					gb.alleles(galleles);
					gb.name(createName.apply(var.sampleName,var.file_index1 -1));
					gb.GQ(var.gq);
					
					gb.attribute(GenpotypeChangedKey,samplesModified.contains(var.sampleName)?1:0);
					gb.attribute(GenpotypeCreated,samplesCreates.contains(var.sampleName)?1:0);
					
					genotypes.add(gb.make());
					}
				b.genotypes(genotypes);
				b.id(first.id);
				
				
				if(!(samplesModified.isEmpty() && samplesCreates.isEmpty()))
					{
					final Set<String> set2 = new TreeSet<String>(samplesModified);
					set2.addAll(samplesCreates);
					b.attribute(GenpotypeDiff, set2.toArray());
					b.filter("DISCORDANCE");
					}
				else
					{
					b.passFilters();
					}
				
				if(!only_print_modified || !(samplesModified.isEmpty() && samplesCreates.isEmpty()))
					{
					w.add(b.make());
					}
						
				}
			iter.close();
			
			w.close();
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(variants!=null) try {variants.cleanup();}catch(Exception err) {}
			}
		return 0;
		}


	public static void main(final String[] args)
		{
		new VCFCompareGT().instanceMainWithExit(args);
		}

}
