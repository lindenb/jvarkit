/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2017 creation

*/
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Example:

content of circos.conf 

```
karyotype = /commun/data/packages/circos/circos-0.69-2/data/karyotype/karyotype.human.hg19.txt
chromosomes_units = 1000000


<<include ideogram.conf>>
<<include ticks.conf>>

<links>
<link>
file          = __FILENAME__
radius        = 0.99r
bezier_radius = 0r
color         = black_a4
thickness     = 2
</link>
</links>


<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>

```

content of ideogram.conf 

```
<ideogram>

<spacing>
## espace entre les contig
default = 0r
</spacing>

radius    = 0.95r
thickness = 30p
fill      = yes

## stroke du radius
stroke_color     = red
stroke_thickness = 3p

# show labels (chromosome names )
show_label       = yes
label_font       = default 
label_radius     = 1r + 75p
label_size       = 30
label_parallel   = yes


</ideogram>
```

content of ticks.conf

```
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p

# the tick label is derived by multiplying the tick position
# by 'multiplier' and casting it in 'format':
#
# sprintf(format,position*multiplier)
#

multiplier       = 1e-6

# %d   - integer
# %f   - float
# %.1f - float with one decimal
# %.2f - float with two decimals
#
# for other formats, see http://perldoc.perl.org/functions/sprintf.html

format           = %d

<tick>
spacing        = 5u
size           = 10p
</tick>

<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>

</ticks>
```

generate circos:


```
$ java -jar dist/samscansplitreads.jar -msr 2 input.bam > input.dat
$ cut -f1-7 input.dat | awk '{OFS="\t";f=$7/50.0;if(f>1.0) f=1.0;if(f<0.5) f=0.5;$1=sprintf("hs%s",$1);$4=sprintf("hs%s",$4);$7=sprintf("thickness=%s,color=%s",f,($1==$4?"blue":"red"));print;}' > tmp.txt
cat circos.conf | sed 's%__FILENAME__%tmp.txt%' >  tmp.txt.conf
circos-0.69-2/bin/circos  -outputdir ./  -outputfile output  -conf  tmp.txt.conf
```



END_DOC
*/
@Program(name="samscansplitreads",description="scan split reads",keywords={"sam","sv","splitreads"})
	public class SamScanSplitReads extends Launcher {
	private static long ID_GENERATOR=0L;
	private static final Logger LOG = Logger.build(SamScanSplitReads.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-x","--extend"},description="extends interval by 'x' pb before merging.")
	private int extentd=20;
	
	@Parameter(names={"-F","--format"},
			description="Output format. if 'vcf', will save the file as a vcf file",
			hidden=true /* this option is not mature */
			)
	private String outputFormat="txt";
	@Parameter(names={"--defaultSampleName"},description="Default Sample name if not read group")
	private String defaultSampleName="UNDEFINED";
	
	@Parameter(names={"--normalize"},description="Optional. Normalize count to this value. e.g '1' ")
	private Integer nomalizeReadCount=null;

	@Parameter(names={"-msr","--minSupportingReads"},description="Minimal number of supporting reads.")
	private int minSupportingReads=0;

	private Map<String,IntervalTreeMap<Set<Arc>>> sample2database = new HashMap<>();
	
	private final Comparator<SAMRecord> coordinateComparator=new Comparator<SAMRecord>()
		{
	    @Override
	    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
	        final String ref1 = samRecord1.getReferenceName();
	        final String ref2 = samRecord2.getReferenceName();
	        final int cmp = ref1.compareTo(ref2);
	        if (cmp != 0) return cmp;
	        return samRecord1.getAlignmentStart() - samRecord2.getAlignmentStart();
	    }
	
		};
		
		
	
	
	private static class Arc
		{
		long id;
		int countSupportingReads=0;
		Interval intervalFrom;
		Interval intervalTo;
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (int) (id ^ (id >>> 32));
			return result;
			}
		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			final Arc other = (Arc) obj;
			return id == other.id;
			}
		boolean intersects(final Arc arc) {
			return intersects(arc.intervalFrom,arc.intervalTo);
			}
		boolean intersects(final Interval i1,final Interval i2) {
			if( this.intervalFrom.intersects(i1) &&
				    this.intervalTo.intersects(i2) ) return true;
			
			if( this.intervalFrom.intersects(i2) &&
				    this.intervalTo.intersects(i1) ) return true;
			return false;
			}
		@Override
		public String toString() {
			return "[" + intervalFrom.toString()+" -> "+intervalTo.toString()+"]";
			}
		}
	
	private Interval extendInterval(final Interval i1) {
		if(this.extentd<=0) return i1;
		return new Interval(i1.getContig(),
				Math.max(i1.getStart()-this.extentd,0),
				i1.getEnd()+this.extentd
				);
		}
	
	private Interval merge(Interval i1,Interval i2) {
		if(!i1.intersects(i2)) throw new IllegalArgumentException(""+i1+"/"+i2);
		return new Interval(i1.getContig(),
				Math.min(i1.getStart(), i2.getStart()),
				Math.max(i1.getEnd(), i2.getEnd())
				);
		}
	
	private void analyseSamPair(final IntervalTreeMap<Set<Arc>> database,final SAMRecord rec1,final SAMRecord rec2)
		{
		int diff = coordinateComparator.compare(rec1, rec2);
		if(diff==0) return;
		if(diff>0) {
			analyseSamPair(database,rec2,rec1);
			return;
			}
		
		final Interval interval1 = extendInterval(new Interval(rec1.getReferenceName(),rec1.getAlignmentStart(), rec1.getAlignmentEnd()));
		final Interval interval2 = extendInterval(new Interval(rec2.getReferenceName(),rec2.getAlignmentStart(), rec2.getAlignmentEnd()));
		if( interval1.intersects(interval2)) return;
		
		
		final Set<Arc> merge=new HashSet<>();
		final Collection<Set<Arc>> col= database.getOverlapping(interval1);
		for(final Set<Arc> arcset:col)
			{
			for(final Arc arc:arcset)
				{
				if(!arc.intersects(interval1,interval2)) continue;
				merge.add(arc);
				}
			}
		
		if(!merge.isEmpty())
			{
			//remove from database
			for(final Set<Arc> arcset:col)
				{
				arcset.removeAll(merge);
				}
			final List<Arc> mergeList=new ArrayList<>(merge);
			//merge new arc with at least one arc in database
			boolean check_one_overlap=false;
			for(final Arc arc: mergeList)
				{
				if(!arc.intervalFrom.intersects(interval1)) continue;
				if(!arc.intervalTo.intersects(interval2)) continue;
				arc.intervalFrom = merge(arc.intervalFrom,interval1);
				arc.intervalTo = merge(arc.intervalTo,interval2);
				arc.countSupportingReads++;
				check_one_overlap=true;
				break;
				}
			if(!check_one_overlap) throw new IllegalStateException();
			int x=0;
			while(x+1<mergeList.size())
				{
				int y=x+1;
				while(y<mergeList.size())
					{
					final Arc arcx = mergeList.get(x);
					final Arc arcy = mergeList.get(y);
					if(arcx.intersects(arcy))
						{
						LOG.debug("mergin "+arcx+" "+arcy+" "+ID_GENERATOR);
						
						arcx.intervalFrom = merge(arcx.intervalFrom,arcy.intervalFrom);
						arcx.intervalTo = merge(arcx.intervalTo,arcy.intervalTo);
						arcx.countSupportingReads+=arcx.countSupportingReads;
						mergeList.remove(y);
						}
					else
						{
						++y;
						}
					}
				
				++x;
				}
			for(final Arc arc: mergeList)
				{
				Set<Arc> arcset = database.get(arc.intervalFrom);
				if(arcset==null) {
					arcset=new HashSet<>();
					 database.put(arc.intervalFrom,arcset);
					}
				arcset.add(arc);
				}
			}
		else
			{
			final Arc arc = new Arc();
			arc.id = ++ID_GENERATOR;
			arc.countSupportingReads=1;
			arc.intervalFrom = interval1;
			arc.intervalTo = interval2;
			Set<Arc> arcset = database.get(interval1);
			if(arcset==null) {
				arcset=new HashSet<>();
				database.put(interval1,arcset);
				}
			arcset.add(arc);
			if(ID_GENERATOR%1000==0) LOG.info(ID_GENERATOR);
			}
		}
	
	private void analyseSamRecord(final SAMRecord rec) {
		if(rec.getReadUnmappedFlag()) return;
		if(rec.getReadFailsVendorQualityCheckFlag()) return;
		if(rec.isSecondaryOrSupplementary()) return;
		if(rec.getDuplicateReadFlag()) return;
		
		final List<SAMRecord> others= SAMUtils.getOtherCanonicalAlignments(rec);
		if(others.isEmpty()) return;
		String sample=this.defaultSampleName;
		final SAMReadGroupRecord g=rec.getReadGroup();
		if(g!=null) {
			final String sa = g.getSample();
			if(sa!=null) sample=sa;
			}
		
		IntervalTreeMap<Set<Arc>> database = this.sample2database.get(sample);
		if(database==null) {
			database=new IntervalTreeMap<Set<Arc>>();
			this.sample2database.put(sample, database);
			}
		
		for(final SAMRecord other:others)
			{
			analyseSamPair(database,rec,other);
			}
		}
	
		private void scanFile(SamReader r) {
			final SAMSequenceDictionaryProgress progess= new SAMSequenceDictionaryProgress(r.getFileHeader());
			final SAMRecordIterator iter = r.iterator();
			while(iter.hasNext())
				{
				analyseSamRecord(progess.watch(iter.next()));
				}
			progess.finish();
			iter.close();
			}
	
		private void saveAsVcf(Set<String> sampleNames,SAMSequenceDictionary dict) throws IOException {
			 final Function<String,Integer> contig2tid = S -> {
				 int i= dict.getSequenceIndex(S);
				 if( i == -1 ) throw new IllegalArgumentException("cannot find contig "+S+" in dictionary");
				 return i;
			 }; 
			
			 final Comparator<Interval> intervalComparator=new Comparator<Interval>()
				{
			    @Override
			    public int compare(final Interval r1, final Interval r2) {
			        final int cmp = contig2tid.apply(r1.getContig())-contig2tid.apply(r2.getContig());
			        if (cmp != 0) return cmp;
			        return r1.getStart() - r2.getStart();
			    	}
				};

			
			final Allele REF = Allele.create("N", true);
			
			final Set<VCFHeaderLine> meta=new HashSet<>();
			VCFStandardHeaderLines.addStandardFormatLines(meta,false,
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.DEPTH_KEY
					);
			sampleNames.addAll(this.sample2database.keySet());
			if(sampleNames.isEmpty()) sampleNames.add(this.defaultSampleName);
			VCFHeader header=new VCFHeader(meta,sampleNames);
			header.setSequenceDictionary(dict);
			VariantContextWriter vcw = super.openVariantContextWriter(outputFile);
			vcw.writeHeader(header);
			final List<Arc> all_arcs = new ArrayList<>();
			
			for(final String sample: this.sample2database.keySet())
				{
				final IntervalTreeMap<Set<Arc>> database = this.sample2database.get(sample);
				for(final Interval interval: database.keySet()) {
					for(final Arc arc: database.get(interval))
						{
						all_arcs.add(arc);
						}
					}
				}
			
			Collections.sort(all_arcs, (A1,A2)->intervalComparator.compare(A1.intervalFrom, A2.intervalFrom));
			
			for(final Arc row:all_arcs)
				{
				final List<Genotype> genotypes = new ArrayList<>();
				final Set<Allele> alleles =  new HashSet<>();
				for(final String sample: this.sample2database.keySet())
					{
					final IntervalTreeMap<Set<Arc>> database = this.sample2database.get(sample);
						for(final Arc arc: database.get(row.intervalFrom))
							{
							if(!arc.equals(row)) continue;
							final Allele alt = Allele.create(
									new StringBuilder().append("<").
									append(arc.intervalFrom.getContig()).append(":").append(arc.intervalFrom.getStart()).append("-").append(arc.intervalFrom.getEnd()).
									append("|").
									append(arc.intervalTo.getContig()).append(":").append(arc.intervalTo.getStart()).append("-").append(arc.intervalTo.getEnd()).
									append(">").toString()
									, false);
							alleles.add(alt);
							final Genotype g = new GenotypeBuilder(sample).alleles(Collections.singletonList(alt)).DP(arc.countSupportingReads).make();
							genotypes.add(g);
							}
					}
				alleles.add(REF);
				VariantContextBuilder vcb=new VariantContextBuilder().
						chr(row.intervalFrom.getContig()).
						start(row.intervalFrom.getStart()).
						stop(row.intervalFrom.getStart()).
						alleles(alleles).
						genotypes(genotypes);
				vcw.add(vcb.make());				
				}
			vcw.close();
			}

		
		private void saveAsText() throws IOException {
			PrintWriter out = super.openFileOrStdoutAsPrintWriter(outputFile);
			
			double maxCount=0.0;
			for(final String sample: this.sample2database.keySet())
				{
				final IntervalTreeMap<Set<Arc>> database = this.sample2database.get(sample);
				for(final Interval interval: database.keySet()) {
					for(final Arc arc: database.get(interval))
						{
						maxCount=Math.max(arc.countSupportingReads, maxCount);
						}
					}
				}
			
			for(final String sample: this.sample2database.keySet())
				{
				final IntervalTreeMap<Set<Arc>> database = this.sample2database.get(sample);
				for(final Interval interval: database.keySet()) {
					for(final Arc arc: database.get(interval))
						{
						if(arc.countSupportingReads<this.minSupportingReads) continue;
						out.print(arc.intervalFrom.getContig());
						out.print("\t");
						out.print(arc.intervalFrom.getStart()-1);
						out.print("\t");
						out.print(arc.intervalFrom.getEnd());
						out.print("\t");
						out.print(arc.intervalTo.getContig());
						out.print("\t");
						out.print(arc.intervalTo.getStart()-1);
						out.print("\t");
						out.print(arc.intervalTo.getEnd());
						out.print("\t");
						if(this.nomalizeReadCount==null) {
							out.print(arc.countSupportingReads);
							}
						else
							{
							out.print(this.nomalizeReadCount.doubleValue()*(arc.countSupportingReads/maxCount));
							}
						out.print("\t");
						out.print(sample);
						out.println();
						}
					}
				}
			out.flush();
			out.close();
		}
		
		private Set<String> samples(SAMFileHeader header )
			{
			return header.getReadGroups().stream().
					map(G->G.getSample()).filter(S->S!=null).
					collect(Collectors.toSet());
			}
		
		@Override
		public int doWork(final List<String> args) {
		 	 SamReader r = null;
		 	 SAMSequenceDictionary dic=null;
		 	 final Set<String> sampleNames=new TreeSet<>();
			try {
				if(args.isEmpty()) {
					LOG.info("read stdin");
					r= openSamReader(null);
					sampleNames.addAll(samples(r.getFileHeader()));
					dic = r.getFileHeader().getSequenceDictionary();
					if(dic==null) {
						LOG.error("SAM input is missing a dictionary");
						return -1;
						}
					scanFile(r);
					r.close();
					r=null;
					}
				else for(final String filename:args)
					{
					LOG.info("read "+filename);
					r= openSamReader(filename);
					sampleNames.addAll(samples(r.getFileHeader()));
					final SAMSequenceDictionary dict2 = r.getFileHeader().getSequenceDictionary();
					if(dict2==null) {
						LOG.error("SAM input is missing a dictionary");
						return -1;
						}
					else if(dic==null)
						{
						dic=dict2;
						}
					else if(!SequenceUtil.areSequenceDictionariesEqual(dic, dict2))
						{
						LOG.error("incompatibles sequences dictionaries");
						return -1;
						}
					scanFile(r);
					r.close();
					r=null;
					}
				
				if("vcf".equalsIgnoreCase(this.outputFormat) || (this.outputFile!=null && (this.outputFile.getName().endsWith(".vcf") || this.outputFile.getName().endsWith(".vcf.gz")))) {
					saveAsVcf(sampleNames,dic);
					}
				else
					{
					saveAsText();
					}
				
				return 0;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);	
			}
		}
	
	public static void main(String[] args) {
		new SamScanSplitReads().instanceMainWithExit(args);
	}
	// 
}
