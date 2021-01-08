/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.nio.file.Path;
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
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
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

with bed files:

```
find . -type f -name "*.bed.gz" > jeter.list
java -jar ${JVARKIT_DIST}/mergecnvnator.jar --input-type bed 
-R ref.fasta jeter.list | bcftools sort -T . -O b -o 20201130.merge.bcf
```


## See also

  * https://github.com/abyzovlab/CNVnator

END_DOC

 */
@Program(name="mergecnvnator",
description="Merge CNVNator results",
keywords= {"cnv","indel","cnvnator"},
biostars={472699},
modificationDate="20201201",
creationDate="20181003"
)
public class MergeCnvNator extends Launcher{
	private static final Logger LOG = Logger.build(MergeCnvNator.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-R","-reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path dictRefFile =  null;
	@Parameter(names={"-r","--ratio"},description="Two intervals are the same if they both have more or equals of this fraction of length in common. " + FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double region_are_same_ratio=0.75;
	@Parameter(names={"--input-type"},description="Input type. Bed type is like cnvnator BUT the 4 first columns are chrom,start,end,sample-name.")
	private InputType inputType = InputType.cnvnator;
	@Parameter(names={"--do-no-reuse"},description="Do not reuse a CNV if it was already used. Undocumented")
	private boolean uniq_use_of_one_cnv =false;
	@Parameter(names={"--one-cnv-type"},description="Only one CNV type (del/dup) per variant.")
	private boolean one_cnv_type =false;
	@Parameter(names={"--treshold"},description="Treshold between HET and HOM genotypes." + FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double treshold =0.2;
	@Parameter(names={"--bed"},description="limit to the calls overlaping that bed FILE")
	private Path includeBed=null;
	@Parameter(names={"--hom-ref"},description="generate HOM_REF instead of NO_CALL for missing genotypes.")
	private boolean hom_ref_instead_of_nocall = false;
	@Parameter(names={"--max-cnv-size"},description="Skip CNVs having a length > 'x'. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int max_cnv_size = 0;
	@Parameter(names={"--ovr"},description="When counting the number of other CNVs overlapping the current CNV in INFO/OV. Just report those having an overlap of at most 'x' with the current variant. The idea is to ignore the big CNVs overlapping the variant." + FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double min_overlapping_ratio=1.0;
	
	
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();
	
	
	private final Allele REF_ALLELE = Allele.create("N", true);
	private final Allele DEL_ALLELE = Allele.create("<DEL>", false);
	private final Allele DUP_ALLELE = Allele.create("<DUP>", false);

	private enum InputType {
		cnvnator,
		bed
	}
	
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
		CNVNatorInterval(final List<String> tokens) {
			this.type = CnvType.valueOf(tokens.get(0));
			final int col = tokens.get(1).indexOf(":");
			if(col<=0) throw new IllegalArgumentException("no colon in "+tokens.get(1));
			this.contig = tokens.get(1).substring(0, col);
			final int hyphen = tokens.get(1).indexOf("-",col+1);
			if(hyphen<0) throw new IllegalArgumentException("no hyphen in "+tokens.get(1));
			this.start= Integer.parseInt(tokens.get(1).substring(col+1, hyphen));
			this.end= Integer.parseInt(tokens.get(1).substring(hyphen+1));
			if(this.start>=this.end)  throw new IllegalArgumentException("bad start/end in "+tokens.get(1));
			}
		CNVNatorInterval(String contig,int start,int end,CnvType type) {
			this.type = type;
			this.contig = contig;
			this.start = start;
			this.end = end;
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
		@SuppressWarnings("unused")
		final Integer pe;
		boolean echoed_flag = false;
		
		CnvNatorCall(String sample /** may be null if input is bed*/, List<String> tokens) {
			final int shift;
			if(sample==null) {
				shift=4;
				this.interval =new CNVNatorInterval(
						tokens.get(0),
						Integer.parseInt(tokens.get(1))+1 /* +1 because BED */,
						Integer.parseInt(tokens.get(2)),
						CnvType.valueOf(tokens.get(4))
						);
				this.sample = tokens.get(3);
				}
			else
				{
				shift = 0;
				this.sample = sample;
				this.interval =new CNVNatorInterval(tokens);
				}
			final Function<Integer,Double> toDbl = IDX->{
				if(IDX>=tokens.size()) return null;
				final String s = tokens.get(IDX);
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
			this.size = toDbl.apply(2 + shift).intValue();//got scientific notation in this column
				
			this.normalized_RD = toDbl.apply(3 + shift);
			this.e_val_1 = toDbl.apply(4 + shift);
			this.e_val_2 = toDbl.apply(5 + shift);
			this.e_val_3 = toDbl.apply(6 + shift);
			this.e_val_4 =toDbl.apply(7 + shift);
			this.q0 = toDbl.apply(8 + shift);
			if(tokens.size()>9+shift ) {
				this.pe = new Integer(tokens.get(9 + shift));
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
		final Interval A = new Interval(a);
		final Interval B = new Interval(b);
		if(!A.intersects(B)) return false;// can happen if baseCall below overlap both
		return testOverlapping2(A,B) &&  testOverlapping2(B,A);
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.region_are_same_ratio<=0 || this.region_are_same_ratio>1) {
			LOG.error("bad region_are_same_ratio :" +this.region_are_same_ratio);
			return -1;
			}
		if(this.treshold< 0 || this.treshold >= 0.5) {
			LOG.error("Bad treshold 0< "+this.treshold+" < 0.5");
			return -1;
			}
		
		try {
			final List<Path> inputs = IOUtils.unrollPaths(args);
			if(inputs.isEmpty()) {
				LOG.error("input is empty");
				return -1;
				}
			final SAMSequenceDictionary dict;
			final ContigNameConverter contigNameConverter;
			if(this.dictRefFile!=null)
				{
				dict = SAMSequenceDictionaryExtractor.extractDictionary(this.dictRefFile);
				contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
				}
			else
				{	
				contigNameConverter = null;
				dict=null;
				}
			
			final IntervalTreeMap<Boolean> limitBedIntervalTreeMap;
			if(this.includeBed!=null) {
				limitBedIntervalTreeMap = new IntervalTreeMap<>();
				final BedLineCodec codec=new BedLineCodec();
				
				try(BufferedReader br= IOUtils.openPathForBufferedReading(this.includeBed)) {
					br.lines().
						filter(S->!StringUtils.isBlank(S) && !BedLine.isBedHeader(S)).
						map(L->codec.decode(L)).
						filter(B->B!=null).
						forEach(B->{
							final String ctg = contigNameConverter==null?B.getContig():contigNameConverter.apply(B.getContig());
							if(StringUtils.isBlank(ctg)) return;
							limitBedIntervalTreeMap.put(new Interval(ctg,B.getStart(),B.getEnd()),Boolean.TRUE);
						});
					}
				}
			else {
				limitBedIntervalTreeMap = null;
				}
			
			final Set<CNVNatorInterval> intervals_set = new HashSet<>();
			final IntervalTreeMap<List<CnvNatorCall>> all_calls = new IntervalTreeMap<>();
			final Set<String> all_samples = new TreeSet<>();
			
			for(final Path input: inputs) {
				final String fileSample;
				LOG.info("Reading "+input);
				if(this.inputType.equals(InputType.cnvnator)) {
					fileSample = IOUtils.getFilenameWithoutCommonSuffixes(input);
					all_samples.add(fileSample);
					}
				else
					{
					fileSample = null;
					}
				final BufferedReader br = IOUtils.openPathForBufferedReading(input);
				String line;
				while((line=br.readLine())!=null) {
					if(StringUtil.isBlank(line) || line.startsWith("#")) continue;
					final String tokens[]= CharSplitter.TAB.split(line);
					final CnvNatorCall call;
					if(this.inputType.equals(InputType.cnvnator)) {
						call = new CnvNatorCall(fileSample, Arrays.asList(tokens));
						}
					else
						{
						call = new CnvNatorCall(null,Arrays.asList(tokens));
						all_samples.add(call.sample);
						}
					if(limitBedIntervalTreeMap!=null && !limitBedIntervalTreeMap.containsOverlapping(call)) {
						continue;
						}
					if(call.size > this.max_cnv_size) {
						continue;
						}
					if(dict!=null && dict.getSequence(call.getContig())==null)
						{
						LOG.warn("skipping "+line+" because contig "+call.getContig()+" is not defined in dictionary");
						continue;
						}
					intervals_set.add(call.interval);
					final Interval key = new Interval(call);
					List<CnvNatorCall> callList = all_calls.get(key);
					if(callList==null) {
						callList = new ArrayList<>();
						all_calls.put(key, callList);
						}
					callList.add(call);
					}
				br.close();
			}
			
			// contig comparator
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
			
			// call comparator
			final Comparator<CNVNatorInterval> comparator = (A,B)->{
				int i = contigComparator.compare(A.getContig(), B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				i = Integer.compare(A.getEnd(), B.getEnd());
				if(i!=0) return i;
				if( one_cnv_type) {
					i = A.type.compareTo(B.type);
					}
				return i;
			};
			final List<CNVNatorInterval> intervals_list = 
					intervals_set.stream().
					sorted(comparator).
					collect(Collectors.toList());
			
			LOG.info("Identified "+ intervals_list.size()+ " intervals");
			
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
			metadata.add(new VCFInfoHeaderLine("NSAMPLES",
					1,
					VCFHeaderLineType.Integer,
					"Number of Samples carrying the CNV"
					));
			
			metadata.add(new VCFInfoHeaderLine(
					"OV",1,
					VCFHeaderLineType.Integer,
					"Number calls overlapping this genotype region, whatever their size ( with at most fraction :"+min_overlapping_ratio+")."
					));

			metadata.add(new VCFFormatHeaderLine(
					"WARN",1,
					VCFHeaderLineType.String,
					"WARNINGS"
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
			JVarkitVersion.getInstance().addMetaData(this, header);
			if(dict!=null) header.setSequenceDictionary(dict);
			
			try(VariantContextWriter out =  this.writingVariants.dictionary(dict).open(this.outputFile)) {
				long id_generator=0L;
				out.writeHeader(header);
				final ProgressFactory.Watcher<CNVNatorInterval> progress = ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
				String prevContig=null;
				for(final CNVNatorInterval interval:intervals_list)
					{
					final Set<String> warnings = new HashSet<>();
					progress.apply(interval);
					
					if(!interval.getContig().equals(prevContig)) {
						//cleanup memory, increase speed ?
						final String finalPrevCtg = prevContig;
						all_calls.keySet().removeIf(R->R.getContig().equals(finalPrevCtg));
						prevContig = interval.getContig();
						}
					
					// find all call overlapping this call
					final List<CnvNatorCall> overlappingList = all_calls.
						getOverlapping(interval).
						stream().
						flatMap(L->L.stream()).
						filter(C->!one_cnv_type || C.interval.type.equals(interval.type)).
						filter(C->!uniq_use_of_one_cnv || !C.echoed_flag).
						collect(Collectors.toList());
					
					if(overlappingList.isEmpty()) continue;
					
					//main call for this interval
					final CnvNatorCall baseCall = overlappingList.stream().
							filter(C->C.interval.equals(interval)).
							findFirst().
							orElse(null);
					
					if(baseCall==null)
						{
						continue;
						}
					// filter calls matching the overlapping with baseCall
					final List<CnvNatorCall> callsToPrint = overlappingList.stream().
							filter(C->testOverlapping(C,baseCall)).
							collect(Collectors.toList())
							;
					
					if(callsToPrint.isEmpty()) {
						throw new IllegalStateException();
						}
					
					
					final VariantContextBuilder vcb= new VariantContextBuilder();
					final Set<Allele> altAlleles = new HashSet<>();
					vcb.chr(baseCall.getContig());
					vcb.start(baseCall.getStart());
					vcb.stop(baseCall.getEnd());
					vcb.attribute(VCFConstants.END_KEY, baseCall.getEnd());
					vcb.attribute("SVLEN", baseCall.size);
					vcb.attribute("IMPRECISE",true);
					
					final Map<String,Genotype> sample2gt = new HashMap<>(callsToPrint.size());
					for(final CnvNatorCall call: callsToPrint)
						{
						if(sample2gt.containsKey(call.sample))
							{
							warnings.add("SAMPLE_"+call.sample+"_MULTIPLE");
							LOG.warn("Sample "+call.sample+" exits twice at the same loc " + call+" could be two small SV overlapping a big one.");
							continue;
							}
						if(uniq_use_of_one_cnv) call.echoed_flag  = true;
						
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
							call.normalized_RD!=null && call.normalized_RD < this.treshold) {
						    gb.alleles(Arrays.asList(DEL_ALLELE,DEL_ALLELE));
						    altAlleles.add(DEL_ALLELE);
						    gb.attribute("CN", 0);
							}
						else if(call.interval.type==CnvType.deletion &&
							call.normalized_RD!=null && call.normalized_RD >= this.treshold)
							{
							gb.alleles(Arrays.asList(REF_ALLELE,DEL_ALLELE));
							altAlleles.add(DEL_ALLELE);
							gb.attribute("CN", 1);
							}
						else if(call.interval.type==CnvType.duplication &&
							call.normalized_RD!=null && call.normalized_RD <= (1.5 + this.treshold))
							{
							gb.alleles(Arrays.asList(REF_ALLELE,DUP_ALLELE));
							altAlleles.add(DUP_ALLELE);
							gb.attribute("CN",2);
							}
						else if(call.interval.type==CnvType.duplication &&
							call.normalized_RD!=null && call.normalized_RD > (1.5 + this.treshold))
							{
							gb.alleles(Arrays.asList(DUP_ALLELE,DUP_ALLELE));
							gb.attribute("CN",9999);
							altAlleles.add(DUP_ALLELE);
							}
						else if(hom_ref_instead_of_nocall) {
							gb.alleles(Arrays.asList(REF_ALLELE,REF_ALLELE));
							}
						else
							{
							gb.alleles(Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
							}
						
						
						sample2gt.put(call.sample,gb.make());
						
						}
					
					
					if(altAlleles.isEmpty()) {
						vcb.attribute(VCFConstants.SVTYPE, "UNDEFINED");
						}
					if(altAlleles.size()==1 && altAlleles.contains(DEL_ALLELE)) {
						vcb.attribute(VCFConstants.SVTYPE, "DEL");
						}
					else if(altAlleles.size()==1 && altAlleles.contains(DUP_ALLELE)) {
						vcb.attribute(VCFConstants.SVTYPE, "DUP");
						}
					else
						{
						vcb.attribute(VCFConstants.SVTYPE, "MIXED");
						}
					final List<Allele> alleles = new ArrayList<>();
					alleles.add(REF_ALLELE);
					alleles.addAll(altAlleles);
					
					vcb.id(String.format("cnv%05d",(++id_generator)));
					vcb.alleles(alleles);
					vcb.attribute("SAMPLES",new ArrayList<>(sample2gt.keySet()));
					vcb.attribute("NSAMPLES",sample2gt.size());
					if(!warnings.isEmpty()) {
						vcb.attribute("WARN",new ArrayList<>(warnings));
						}
					
					vcb.attribute(
							"OV",
							 all_calls.
								getOverlapping(baseCall).
								stream().
								flatMap(L->L.stream()).
								filter(L->{
									final Interval r1 = new Interval(L);
									final Interval r2 = new Interval(baseCall);
									if(!r1.overlaps(r2)) return false;
									final double L1 = r1.getIntersectionLength(r2);
									if(  (L1 /r1.length() )  < min_overlapping_ratio) return false;
									return true;
									}).
								count()
							);
					
					vcb.genotypes(sample2gt.values());
					out.add(vcb.make());
					}
				progress.close();
				}
			
			if(uniq_use_of_one_cnv) {
				all_calls.values().stream().
					flatMap(C->C.stream()).
					filter(C->!C.echoed_flag).
					forEach(C->LOG.warn("Bug: Not printed "+C));
				}
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
		}
		}
	
	
public static void main(final String[] args) {
	new MergeCnvNator().instanceMainWithExit(args);
}
}
