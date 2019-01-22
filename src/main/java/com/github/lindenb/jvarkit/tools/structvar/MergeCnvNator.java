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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
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
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## input

Input is the tabular output of CNVNator

```
(...)
deletion     chr2:1-2649000        2.649e+06  0        6.01633e-14  0            6.02087e-14  0           -1
duplication  chr2:3712001-3721000  9000       1.89036  0.0274362    5.17568e-42  0.137369     3.6838e-64  0.00821444
(...)
```

The name of each sample is the `basename` of the file, before the first `.`

a list of paths can be specified if the only input file ends with '.list' 


## Example

```
find DIR1 DIR2 -type f -name "*.tsv" > in.list
java -jar dist/mergecnvnator.jar -R ref.fasta in.list > out.vcf

(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample7	Sample8	Sample9
chr1	1	.	N	<DEL>	.	.	END=10000;IMPRECISE;SVLEN=10000;SVTYPE=DEL	GT:CN:P1:P2:P3:P4:Q0:RD	1/1:0:1.594e-11:0.00:1.992e-11:0.00:-1.000e+00:0.00	1/1:0:1.594e-11:0.00:1.992e-11:0.00:-1.000e+00:0.00	1/1:0:1.594e-11:0.00:1.992e-11:0.00:-1.000e+00:0.00
chr1	40001	.	N	<DEL>	.	.	END=48000;IMPRECISE;SVLEN=8000;SVTYPE=DEL	GT:CN:P1:P2:P3:P4:Q0:RD	./.	./.	0/1:1:34.95:5.850e-06:323.37:2.143e-03:0.889:0.603
(...)
```

## See also

  * https://github.com/abyzovlab/CNVnator

END_DOC

 */
@Program(name="mergecnvnator",
description="Merge CNVNator results",
keywords= {"cnv","indel","cnvnator"}
)
public class MergeCnvNator extends Launcher{
	private static final Logger LOG = Logger.build(MergeCnvNator.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-R","-reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File dictRefFile =  null;
	@Parameter(names={"-r","--ratio"},description="two intervals are the same if they both have more or equals of this fraction of length in common")
	private double region_are_same_ratio=0.75;
	
	private final Allele REF_ALLELE = Allele.create("N", true);
	private final Allele DEL_ALLELE = Allele.create("<DEL>", false);
	private final Allele DUP_ALLELE = Allele.create("<DUP>", false);

	private enum CnvType {
		deletion,
		duplication
	} ;
	
	private class CNVNatorInterval implements Locatable
		{
		final CnvType type;
		final String contig;
		final int start;
		final int end;
		CNVNatorInterval(final String tokens[]) {
			this.type = CnvType.valueOf(tokens[0]);
			final int col = tokens[1].indexOf(":");
			if(col<=0) throw new IllegalArgumentException("no colon in "+tokens[1]);
			this.contig = tokens[1].substring(0, col);
			final int hyphen = tokens[1].indexOf("-",col+1);
			if(hyphen<0) throw new IllegalArgumentException("no hyphen in "+tokens[1]);
			this.start= Integer.parseInt(tokens[1].substring(col+1, hyphen));
			this.end= Integer.parseInt(tokens[1].substring(hyphen+1));
			if(this.start>=this.end)  throw new IllegalArgumentException("bad start/end in "+tokens[1]);
			}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + contig.hashCode();
			result = prime * result + end;
			result = prime * result + start;
			result = prime * result + type.hashCode();
			return result;
		}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj)  return true;
			if (obj == null)  return false;
			if (!(obj instanceof CNVNatorInterval)) {
				return false;
				}
			final CNVNatorInterval other = (CNVNatorInterval) obj;
			if (end != other.end) {
				return false;
				}
			if (start != other.start) {
				return false;
				}
			if (type != other.type) {
				return false;
				}
			if (!contig.equals(other.contig)) {
				return false;
				}
			return true;
			}

		@Override
		public String getContig() {
			return this.contig;
			}
		@Override
		public int getStart() {
			return this.start;
			}
		@Override
		public int getEnd() {
			return this.end;
			}
		@Override
		public String toString() {
			return type+"/"+this.contig+":"+this.start+"-"+this.end;
			}
		}
	
	// CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0
	private class CnvNatorCall implements Locatable
		{
		final String sample;
		final CNVNatorInterval interval;
		final int size;
		final Double normalized_RD ;
		final Double e_val_1 ;
		final Double e_val_2 ;
		final Double e_val_3 ;
		final Double e_val_4 ;
		final Double q0;
		final Integer pe;
		boolean echoed_flag = false;
		
		CnvNatorCall(String sample,final String tokens[]) {
			this.sample = sample;
			this.interval =new CNVNatorInterval(tokens);
		
			final Function<Integer,Double> toDbl = IDX->{
				if(IDX>=tokens.length) return null;
				final String s = tokens[IDX];
				if(s.isEmpty()) return null;
				if(s.equals("nan") || s.equals("-nan")) return null;
				try {
					return Double.valueOf(s);
					}
				catch(final NumberFormatException err) {
					LOG.warn(err);
					return null;
					}
				};
			this.size = toDbl.apply(2).intValue();//got scientific notation in this column
				
			this.normalized_RD = toDbl.apply(3);
			this.e_val_1 = toDbl.apply(4);
			this.e_val_2 = toDbl.apply(5);
			this.e_val_3 = toDbl.apply(6);
			this.e_val_4 =toDbl.apply(7);
			this.q0 = toDbl.apply(8);
			if(tokens.length>9 ) {
				this.pe = new Integer(tokens[9]);
				}
			else
				{
				this.pe = null;
				}
			}
		@Override
		public String getContig() {
			return interval.getContig();
			}
		@Override
		public int getStart() {
			return interval.getStart();
			}
		@Override
		public int getEnd() {
			return interval.getEnd();
			}
		@Override
		public String toString() {
			return sample+" "+this.interval;
			}
		}
	
	private boolean testOverlapping2(final Interval A,final Interval B )
		{
		double lenA = A.length();
		double len2 = A.getIntersectionLength(B);
		return len2/lenA >= this.region_are_same_ratio;
		}
	
	private boolean testOverlapping(final Locatable a,final Locatable b )
		{
		final Interval A = new Interval(a.getContig(), a.getStart(), a.getEnd());
		final Interval B = new Interval(b.getContig(), b.getStart(), b.getEnd());
		if(!A.intersects(B)) return false;// can happen if baseCall below overlap both
		return testOverlapping2(A,B) &&  testOverlapping2(B,A);
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.region_are_same_ratio<=0 || this.region_are_same_ratio>1) {
			LOG.error("bad region_are_same_ratio :" +this.region_are_same_ratio);
			return -1;
			}
		VariantContextWriter out = null;
		try {
			final List<File> inputs = IOUtils.unrollFiles2018(args);
			if(inputs.isEmpty()) {
				LOG.error("input is empty");
				return -1;
				}
			final SAMSequenceDictionary dict;
			if(this.dictRefFile!=null)
				{
				dict = SAMSequenceDictionaryExtractor.extractDictionary(this.dictRefFile);
				}
			else
				{	
				dict=null;
				}
			
			final Set<CNVNatorInterval> intervals_set = new HashSet<>();
			final IntervalTreeMap<List<CnvNatorCall>> all_calls = new IntervalTreeMap<>();
			final Set<String> all_samples = new TreeSet<>();
			
			for(final File input: inputs) {
				String sample = input.getName();
				int dot = sample.indexOf('.');
				if(dot>0) sample=sample.substring(0,dot);
				all_samples.add(sample);
				final BufferedReader br = IOUtils.openFileForBufferedReading(input);
				String line;
				while((line=br.readLine())!=null) {
					if(StringUtil.isBlank(line)) continue;
					final String tokens[]= CharSplitter.TAB.split(line);
					final CnvNatorCall call = new CnvNatorCall(sample, tokens);
					if(dict!=null && dict.getSequence(call.getContig())==null)
						{
						LOG.warn("skipping "+line+" because contig "+call.getContig()+" is not defined in dictionary");
						continue;
						}
					intervals_set.add(call.interval);
					final Interval key = new Interval(call.getContig(), call.getStart(), call.getEnd());
					List<CnvNatorCall> callList = all_calls.get(key);
					if(callList==null) {
						callList = new ArrayList<>();
						all_calls.put(key, callList);
						}
					callList.add(call);
				}
				br.close();
			}
			
			final Comparator<String> contigComparator;
			if(dict!=null)
				{
				contigComparator = new ContigDictComparator(dict);
				}
			else
				{
				final SmartComparator smartComparator = new SmartComparator();
				contigComparator = (A,B)-> smartComparator.compare(A, B);
				}
			
			final Comparator<CNVNatorInterval> comparator = (A,B)->{
				int i = contigComparator.compare(A.getContig(), B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				i = Integer.compare(A.getEnd(), B.getEnd());
				if(i!=0) return i;
				i = A.type.compareTo(B.type);
				return i;
			};
			final List<CNVNatorInterval> intervals_list = 
					intervals_set.stream().
					sorted(comparator).
					collect(Collectors.toList());
			
			
			final Set<VCFHeaderLine> metadata = new HashSet<>();
			
			VCFStandardHeaderLines.addStandardFormatLines(metadata, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.GENOTYPE_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(metadata, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.END_KEY
					);
			metadata.add(new VCFInfoHeaderLine("SVLEN",1, VCFHeaderLineType.Integer, "Difference in length between REF and ALT alleles"));
			metadata.add(new VCFInfoHeaderLine("IMPRECISE",0, VCFHeaderLineType.Flag, "Imprecise structural variation"));
			metadata.add(new VCFInfoHeaderLine(VCFConstants.SVTYPE,1, VCFHeaderLineType.String, "Structural variation type"));
			metadata.add(new VCFInfoHeaderLine("SAMPLES",
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Samples carrying the CNV"
					));

			
			metadata.add(new VCFFormatHeaderLine(
					"RD",1,
					VCFHeaderLineType.Float,
					"Normalized RD"
					));
			metadata.add(new VCFFormatHeaderLine(
					"P1",1,
					VCFHeaderLineType.Float,
					"e-val by t-test"
					));
			metadata.add(new VCFFormatHeaderLine(
					"P2",1,
					VCFHeaderLineType.Float,
					"e-val by Gaussian tail"
					));
			metadata.add(new VCFFormatHeaderLine(
					"P3",1,
					VCFHeaderLineType.Float,
					"e-val by t-test (middle)"
					));
			
			metadata.add(new VCFFormatHeaderLine(
					"P4",1,
					VCFHeaderLineType.Float,
					"e-val by Gaussian tail (middle)"
					));
			metadata.add(new VCFFormatHeaderLine(
					"Q0",1,
					VCFHeaderLineType.Float,
					"Fraction of reads with 0 mapping quality"
					));
			
			metadata.add(new VCFFormatHeaderLine(
					"CN",1,
					VCFHeaderLineType.Integer,
					"Copy number genotype for imprecise events"
					));
			metadata.add(new VCFFormatHeaderLine(
					"PE",1,
					VCFHeaderLineType.Integer,
					"Number of paired-ends that support the event"
					));
			
			metadata.add(new VCFFormatHeaderLine(
					"OV",1,
					VCFHeaderLineType.Integer,
					"Number calls (with different sample) overlapping this genotype"
					));
			
			final VCFHeader header = new VCFHeader(
					metadata,
					all_samples
					);
			if(dict!=null) header.setSequenceDictionary(dict);
			
			out =  super.openVariantContextWriter(this.outputFile);
			out.writeHeader(header);
			for(final CNVNatorInterval interval:intervals_list)
				{
				final List<CnvNatorCall> overlappingList = all_calls.
					getOverlapping(interval).
					stream().
					flatMap(L->L.stream()).
					filter(C->C.interval.type.equals(interval.type)).
					filter(C->!C.echoed_flag).
					collect(Collectors.toList());
				
				if(overlappingList.isEmpty()) continue;
				final CnvNatorCall baseCall = overlappingList.stream().
						filter(C->C.interval.equals(interval)).
						findFirst().
						orElse(null);
				if(baseCall==null)
					{
					continue;
					}
				final List<CnvNatorCall> callsToPrint = overlappingList.stream().
						filter(C->testOverlapping(C,baseCall)).
						collect(Collectors.toList())
						;
				
				if(callsToPrint.isEmpty()) throw new IllegalStateException();
				
				final VariantContextBuilder vcb= new VariantContextBuilder();
				vcb.chr(baseCall.getContig());
				vcb.start(baseCall.getStart());
				vcb.stop(baseCall.getEnd());
				vcb.attribute(VCFConstants.END_KEY, baseCall.getEnd());
				vcb.attribute("SVLEN", baseCall.size);
				vcb.attribute("IMPRECISE",true);
				
				switch(baseCall.interval.type)
					{
					case deletion:
						vcb.alleles(Arrays.asList(REF_ALLELE,DEL_ALLELE));
						vcb.attribute(VCFConstants.SVTYPE, "DEL");
						break;
					case duplication:
						vcb.alleles(Arrays.asList(REF_ALLELE,DUP_ALLELE));
						vcb.attribute(VCFConstants.SVTYPE, "DUP");
						break;
					default: throw new IllegalStateException();
					}
				
				
				final Map<String,Genotype> sample2gt = new HashMap<>(callsToPrint.size());
				for(final CnvNatorCall call: callsToPrint)
					{
					if(sample2gt.containsKey(call.sample))
						{
						LOG.warn("Sample "+call.sample+" exits twice at the same loc " + call+" could be two small SV overlapping a big one.");
						continue;
						}
					call.echoed_flag  = true;
					
					final GenotypeBuilder gb = new GenotypeBuilder(call.sample);
					
					if( call.normalized_RD!=null) gb.attribute("RD", call.normalized_RD);
					if( call.e_val_1!=null) gb.attribute("P1", call.e_val_1);
					if( call.e_val_2!=null) gb.attribute("P2", call.e_val_2);
					if( call.e_val_3!=null) gb.attribute("P3", call.e_val_3);
					if( call.e_val_4!=null) gb.attribute("P4", call.e_val_4);
					if( call.q0!=null) gb.attribute("Q0", call.q0);
					gb.attribute("OV",
							all_calls.getOverlapping(call).
								stream().
								flatMap(col->col.stream()).
								filter(C->!C.sample.equals(call.sample)).
								count()
							);
					
					if (call.interval.type==CnvType.deletion &&
						call.normalized_RD!=null && call.normalized_RD <0.20) {
					    gb.alleles(Arrays.asList(DEL_ALLELE,DEL_ALLELE));
					    gb.attribute("CN", 0);
						}
					else if(call.interval.type==CnvType.deletion &&
						call.normalized_RD!=null && call.normalized_RD >=0.20)
						{
						gb.alleles(Arrays.asList(REF_ALLELE,DEL_ALLELE));
						gb.attribute("CN", 1);
						}
					else if(call.interval.type==CnvType.duplication &&
						call.normalized_RD!=null && call.normalized_RD <=1.7)
						{
						gb.alleles(Arrays.asList(REF_ALLELE,DUP_ALLELE));
						gb.attribute("CN",2);
						}
					else if(call.interval.type==CnvType.duplication &&
						call.normalized_RD!=null && call.normalized_RD >1.7)
						{
						gb.alleles(Arrays.asList(DUP_ALLELE,DUP_ALLELE));
						gb.attribute("CN",9999);
						}
					else
						{
						gb.alleles(Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
						}
					
					
					sample2gt.put(call.sample,gb.make());
					
					}
				vcb.attribute("SAMPLES",new ArrayList<>(sample2gt.keySet()));
				vcb.genotypes(sample2gt.values());
				out.add(vcb.make());
				}
			
			all_calls.values().stream().
				flatMap(C->C.stream()).
				filter(C->!C.echoed_flag).
				forEach(C->LOG.warn("Bug: Not printed "+C));
			
			out.close();
			out=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(out);
		}
		}
	
	
public static void main(final String[] args) {
	new MergeCnvNator().instanceMainWithExit(args);
}
}
