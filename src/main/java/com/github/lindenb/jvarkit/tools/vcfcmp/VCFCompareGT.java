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
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfIteratorImpl;

/**
BEGIN_DOC

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

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	private class Variant
		{
		String chrom="";
		String id=VCFConstants.EMPTY_ID_FIELD;
		String ref=VCFConstants.EMPTY_ALLELE;
		int start=-1;
		int end=-1;
		int file_index=0;
		String sampleName="";
		String a1=VCFConstants.EMPTY_ALLELE;
		String a2=VCFConstants.EMPTY_ALLELE;
		int dp=-1;
		int gq=-1;
		}
	
	private class PosComparator implements Comparator<Variant>
		{
		@Override
		public int compare(Variant v1, Variant v2) {
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
		public int compare(Variant v1, Variant v2) {
			int i=super.compare(v1, v2);
			if(i!=0) return i;
			i=v1.sampleName.compareTo(v2.sampleName);
			if(i!=0) return i;
			i=v1.file_index - v2.file_index;
			if(i!=0) return i;
			return 0;
			}
		}
	
	private class VariantCodec extends AbstractDataCodec<Variant>
		{
		@Override
		public Variant decode(DataInputStream dis) throws IOException
			{
			Variant v=new Variant();
			try {
				v.chrom=dis.readUTF();
			} catch (Exception e) {
				return null;
				}
			v.start=dis.readInt();
			v.end=dis.readInt();
			v.id=dis.readUTF();
			v.ref=dis.readUTF();
			v.sampleName=dis.readUTF();
			v.a1=readString(dis);
			v.a2=readString(dis);
			v.file_index=dis.readInt();
			v.dp=dis.readInt();
			v.gq=dis.readInt();
			return v;
			}
		@Override
		public void encode(DataOutputStream dos, Variant v)
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
			dos.writeInt(v.file_index);
			dos.writeInt(v.dp);
			dos.writeInt(v.gq);
			}
		@Override
		public AbstractDataCodec<Variant> clone() {
			return new VariantCodec();
			}
		}
	
	
	private VCFCompareGT()
		{
		}
	


	@Override
	public int doWork(final List<String> args) {
		if(args.isEmpty())
			{
			System.err.println("VCF missing.");
			return -1;
			}
		VariantComparator varcmp=new VariantComparator();
		SortingCollection<Variant> variants = null;
		Set<String> sampleNames=new LinkedHashSet<String>();
		try
			{
			
			variants=SortingCollection.newInstance(
					Variant.class,
					new VariantCodec(),
					varcmp,
					writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpDirectories()
					);
			variants.setDestructiveIteration(true);
			
			Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
			metaData.add(new VCFHeaderLine(getClass().getSimpleName(),"version:"+getVersion()+" command:"+getProgramCommandLine()));
			
			
			for(int i=0;i< args.size();++i)
				{
				File vcfFile=new File(args.get(i));
				LOG.info("Opening "+vcfFile);
				VcfIterator iter= new VcfIteratorImpl(IOUtils.openFileForReading(vcfFile));
				VCFHeader header=iter.getHeader();
				sampleNames.addAll(header.getSampleNamesInOrder());
				
				metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"_"+((i)+1),"File: "+vcfFile.getPath()));

				
				long nLines=0;
				while(iter.hasNext())
					{
					VariantContext var=iter.next();
					
					if(nLines++%10000==0)
						{
						LOG.info(args.get(i)+" "+nLines);
						}
					if(var.getReference()==null || var.getReference().isSymbolic()) continue;
					if(!var.hasGenotypes()) continue;
					for(Genotype genotype:var.getGenotypes())
						{
						Variant rec=new Variant();
						if(!genotype.isAvailable()) continue;
						if(!genotype.isCalled()) continue;
						if(genotype.isNoCall()) continue;
						
						rec.file_index=i+1;
						rec.sampleName=genotype.getSampleName();
						rec.chrom=var.getContig();
						rec.start=var.getStart();
						rec.end=var.getEnd();
						
						rec.ref=var.getReference().getBaseString();
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
						List<Allele> alleles=genotype.getAlleles();
						if(alleles==null) continue;
						if(alleles.size()==1)
							{
							rec.a1=alleles.get(0).getBaseString().toUpperCase();
							rec.a2=rec.a1;
							}
						else if(alleles.size()==2)
							{
							rec.a1=alleles.get(0).getBaseString().toUpperCase();
							rec.a2=alleles.get(1).getBaseString().toUpperCase();
							if(rec.a1.compareTo(rec.a2)>0)
								{
								String tmp=rec.a2;
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
				}
			variants.doneAdding();
		
			LOG.info("Done Adding");
			
			Set<String> newSampleNames=new HashSet<String>();
			for(int i=0;i< args.size();++i)
				{
				for(String sample:sampleNames)
					{
					newSampleNames.add(sample+"_"+((i)+1));
					}
				}
			final String GenpotypeChangedKey="GCH";
			final String GenpotypeCreated="GNW";
			final String GenpotypeDiff="GDF";
			
			metaData.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
			metaData.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Depth"));
			metaData.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Integer, "Qual"));
			metaData.add(new VCFFormatHeaderLine(GenpotypeChangedKey,1,VCFHeaderLineType.Integer, "Changed Genotype"));
			metaData.add(new VCFFormatHeaderLine(GenpotypeCreated,1,VCFHeaderLineType.Integer, "Genotype Created/Deleted"));
			metaData.add(new VCFInfoHeaderLine(GenpotypeDiff,VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String, "Samples with Genotype Difference"));

			
			
			VCFHeader header=new VCFHeader(
					metaData,
					new ArrayList<String>(newSampleNames));
			
			VariantContextWriter w= VCFUtils.createVariantContextWriterToStdout();
			w.writeHeader(header);
			List<Variant> row=new ArrayList<Variant>();
			final PosComparator posCompare=new PosComparator();
			CloseableIterator<Variant> iter=variants.iterator();
			for(;;)
				{
				Variant rec=null;
				if(iter.hasNext())
					{
					rec=iter.next();
					}
				if(rec==null || (!row.isEmpty() && posCompare.compare(row.get(0),rec)!=0))
					{
					if(!row.isEmpty())
						{
						Set<String> samplesModified=new TreeSet<String>();
						Set<String> samplesCreates=new TreeSet<String>();
						Counter<String> samplesSeen=new Counter<String>();
						for(int x=0;x< row.size();++x)
							{
							Variant var1=row.get(x);
							samplesSeen.incr(var1.sampleName);
							for(int y=x+1;y< row.size();++y)
								{
								Variant var2=row.get(y);
								if(!var2.sampleName.equals(var1.sampleName)) continue;
								if(var1.a1.equals(var2.a1) && var1.a2.equals(var2.a2) ) continue;
								samplesModified.add(var1.sampleName);
								}
							}
						
						for(String sampleName:samplesSeen.keySet())
							{
							if(samplesSeen.count(sampleName)!=args.size())
								{
								samplesCreates.add(sampleName);
								}
							}
						
						
						Variant first=row.get(0);
						Set<Allele> alleles=new HashSet<Allele>();
						alleles.add(Allele.create(first.ref, true));
						for(Variant var:row)
							{
							alleles.add(Allele.create(var.a1, var.a1.equalsIgnoreCase(var.ref)));
							alleles.add(Allele.create(var.a2, var.a2.equalsIgnoreCase(var.ref)));
							}
						
						
						VariantContextBuilder b=new VariantContextBuilder(
								getClass().getName(),
								first.chrom,
								first.start,
								first.end,
								alleles
								);
						//build genotypes
						List<Genotype> genotypes=new ArrayList<Genotype>();
						for(Variant var:row)
							{
							//alleles for this genotype
							List<Allele> galleles=new ArrayList<Allele>();
							galleles.add(Allele.create(var.a1, var.a1.equalsIgnoreCase(var.ref)));
							galleles.add(Allele.create(var.a2, var.a2.equalsIgnoreCase(var.ref)));
							
							GenotypeBuilder gb=new GenotypeBuilder();
							gb.DP(var.dp);
							gb.alleles(galleles);
							gb.name(var.sampleName+"_"+var.file_index);
							gb.GQ(var.gq);
							
							gb.attribute(GenpotypeChangedKey,samplesModified.contains(var.sampleName)?1:0);
							gb.attribute(GenpotypeCreated,samplesCreates.contains(var.sampleName)?1:0);
							
							genotypes.add(gb.make());
							
							
							
							}
						b.genotypes(genotypes);
						b.id(first.id);
						
						
						if(!(samplesModified.isEmpty() && samplesCreates.isEmpty()))
							{
							Set<String> set2=new TreeSet<String>(samplesModified);
							set2.addAll(samplesCreates);
							b.attribute(GenpotypeDiff, set2.toArray());
							}
						
						
						if(!only_print_modified || !(samplesModified.isEmpty() && samplesCreates.isEmpty()))
							{
							w.add(b.make());
							}
						row.clear();
						}
					if(rec==null) break;
					}
				row.add(rec);
				}
			iter.close();
			
			w.close();
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(variants!=null) variants.cleanup();
			}
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VCFCompareGT().instanceMainWithExit(args);
		}

}
