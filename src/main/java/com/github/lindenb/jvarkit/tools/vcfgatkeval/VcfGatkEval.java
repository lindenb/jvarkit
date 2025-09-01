/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfgatkeval;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC
 
# Run on one vcf

```
$ bcftools view in.bcf | java -jar dist/jvarkit.jar  vcfgatkeval -o "out1"
```

list output

```
$ ls out1.*
out1.output.filters.txt
out1.output.R 
out1.output.table.txt
```

plot barplots:

```
R --vanilla --no-save < out1.output.R 
```


filters for gatk that can be used using gatk `--arguments_file`

```
cat  out1.output.filters.txt
```

```
-filter
"vc.isSNP() && FS > 21.0"
--filter-name
FS_HIGH_SNP
-filter
"vc.isSNP() && MQ < 60.0"
--filter-name
MQ_LOW_SNP
-filter
"vc.isSNP() && MQRankSum < 0.0"
--filter-name
MQRankSum_LOW_SNP
-filter
"vc.isSNP() && MQRankSum > 0.0"
--filter-name
MQRankSum_HIGH_SNP
-filter
"vc.isSNP() && QD < 1.0"
--filter-name
QD_LOW_SNP
-filter
"vc.isSNP() && ReadPosRankSum < -2.2"
--filter-name
ReadPosRankSum_LOW_SNP
-filter
"vc.isSNP() && ReadPosRankSum > 2.4"
--filter-name
ReadPosRankSum_HIGH_SNP
-filter
"vc.isSNP() && SOR > 3.5"
--filter-name
SOR_HIGH_SNP
```

use with gatk variantFilteration:

```
gatk VariantFiltration -V in.vcf.gz -R reference.fasta -O out.vcf.gz --arguments_file out1.output.filters.txt

```

# Parallelisation

run in parallel
```
$ bcftools view in.bcf chr1 | java -jar dist/jvarkit.jar  vcfgatkeval -o "out1"
$ bcftools view in.bcf chr2 | java -jar dist/jvarkit.jar  vcfgatkeval -o "out2"
```

and then concat:

```
cat out1.output.table.txt out2.output.table.txt  |  java -jar dist/jvarkit.jar  vcfgatkeval --input-type table  -o "out3"
```


END_DOC
*/

@Program(
		name="vcfgatkeval",
		description="Eval/Plot gatk INFO tags for filtering",
		jvarkit_amalgamion = true,
		keywords={"vcf","gatk"},
		creationDate ="20230424",
		modificationDate="20240321"
		)
public class VcfGatkEval extends Launcher {
	private static Logger LOG=Logger.of(VcfGatkEval.class);
	private enum INPUT_TYPE{vcf,table};
	
	@Parameter(names={"-o","--output"},description="filename prefix")
	private String basename ="gatk.eval";
	@Parameter(names={"-I","--input-type"},description="input type. vcf: vcf file or stdin. table: stdin or one or more output of *.output.table.txt")
	private INPUT_TYPE inputType = INPUT_TYPE.vcf;
	@Parameter(names={"-p","--percentile"},description="GATK Filters should be applied to this percentile: f < x < (1.0 -f)")
	private double percentile= 0.025;
	@Parameter(names={"--depth","--with-depth"},description="include INFO/" + VCFConstants.DEPTH_KEY)
	private boolean include_depth = false;

	
	private static class Range {
		final double lowerBound;
		final double upperBound;
		long count_ALL = 0L;
		long count_PASS = 0L;
		Range(final double m,final double M) {
			this.lowerBound = m;
			this.upperBound = M;
			}
		boolean contains(final double f) {
			return this.lowerBound <= f && f < this.upperBound;
			}
		@Override
		public String toString() {
			return "("+lowerBound+"/"+upperBound+"(:PASS~"+count_PASS+"/ALL~"+count_ALL;
			}
		}
		
	private abstract class Annotation {
		final VariantContext.Type variantType;
		final List<Range> ranges = new ArrayList<>();
		protected VCFInfoHeaderLine infoHeaderLine = null;
		private final DecimalFormat decimalFormat;
		Annotation(final VariantContext.Type variantType) {
			this.variantType = variantType;
			String str="#.";
			int p = getPrecision();
			while(p>1) {
				str+="#";
				p=p/10;
				}
			this.decimalFormat = new DecimalFormat(str);
			this.decimalFormat.setRoundingMode(java.math.RoundingMode.CEILING);
			}
		
		abstract String getTag();
		void init(final VCFHeader header) {
			this.infoHeaderLine = header.getInfoHeaderLine(getTag());
			if(this.infoHeaderLine==null) {
				LOG.warning("The following tag was not found in header:"+getTag());
				}
			}
		

		abstract int getPrecision();
		abstract String getSub();

		boolean logX() { return false;}
		

		double round(final double v) {
			return Double.parseDouble(this.decimalFormat.format(v));
			}		

		
		void visit(final VariantContext ctx) {
			if(this.infoHeaderLine==null) return;
			if(!ctx.getType().equals(this.variantType)) return;
			if(!ctx.hasAttribute(getTag())) return;
			try {
			for(String s : ctx.getAttributeAsStringList(getTag(),".") ) {
				if(s.equals(".")) continue;
				final double v0 = Double.parseDouble(s);
				if(Double.isNaN(v0)) continue;
				if(Double.isInfinite(v0)) continue;
				final double v = round(v0); 
				int i=0;
				for(i=0;i< this.ranges.size();i++) {
					final Range r = ranges.get(i);
					if(r.contains(v)) {
						r.count_ALL++;
						if(ctx.isNotFiltered()) {
							r.count_PASS++;
							}
						break;
						}
					}
				if(i!=this.ranges.size()) continue;
				
				final int precision = getPrecision();
				final double v1 = v;
				final double v2 = v1 + 1.0/precision;
				final Range r = new Range(v1,v2);
				if(!r.contains(v)) throw new IllegalStateException(""+r+" "+v);
				r.count_ALL = 1L;
				if(ctx.isNotFiltered()) {
					r.count_PASS=1L;
					}

				if(this.ranges.stream().anyMatch(R->R.lowerBound==r.lowerBound)) {
					System.err.println("r.size "+this.ranges.size());
					for(Range r2: this.ranges) {
						System.err.println("# "+r2+"="+r2.contains(v)+" v="+v);
						}
					throw new IllegalStateException(""+this.ranges+" -- and -- "+r+" i="+i+" v="+v+" v1="+v1+" v2="+v2+" precision="+precision);
					}
		
				this.ranges.add(r);
				Collections.sort(this.ranges,(A,B)->Double.compare(A.lowerBound, B.lowerBound));
				}
			} catch(Throwable err) {
				LOG.severe("Error "+err.getMessage()+" "+ctx+" "+toString());
				throw err;
				}
			}
		
		void saveFilter(PrintWriter out) throws IOException {
			if(this.ranges.isEmpty()) return ;
			final String testType;
			switch(this.variantType) {
				case SNP: testType = "vc.isSNP()&&"; break;//no blank or gatk fails
				case INDEL: testType = "vc.isIndel()&&"; break;//no blank or gatk fails
				default: testType = ""; break;
				}
			final String hasAttribute = "vc.hasAttribute(\""+getTag()+"\")&&";
			OptionalDouble limit = getLowPercentile();
			if(limit.isPresent()) {
				out.println("--filter-expression");
				out.println(testType + hasAttribute +  getTag() +"<"+ limit.getAsDouble()); //no blank or gatk fails
				out.println("--filter-name");
				out.println(getTag()+"_LOW_"+variantType.name());
				}
			limit = getHighPercentile();
			if(limit.isPresent()) {
				out.println("--filter-expression");
				out.println(testType + hasAttribute + getTag() +">"+ limit.getAsDouble());//no blank or gatk fails
				out.println("--filter-name");
				out.println(getTag()+"_HIGH_"+variantType.name());
				}

			}
		
		void saveTable(PrintWriter out) throws IOException {
			out.println(getTag());
			out.println(this.variantType.name());
			out.println(this.ranges.size());
			for(Range r: ranges) {
				out.print(r.lowerBound);
				out.print("\t");
				out.print(r.upperBound);
				out.print("\t");
				out.print(r.count_ALL);
				out.print("\t");
				out.println(r.count_PASS);
				}
			}
		
		void loadTable(BufferedReader br) throws IOException {
			String line = br.readLine();
			if(line==null) throw new IOException("count missing for "+getTag());
			int n = Integer.parseInt(line);
			for(int i=0;i< n;i++) {
				line = br.readLine();
				if(line==null) throw new IOException("line["+i+"] missing for "+getTag());
				final String[] tokens = CharSplitter.TAB.split(line);
				final double m = round(Double.parseDouble(tokens[0]));
				final double M = round(Double.parseDouble(tokens[1]));
				final long count_ALL = Long.parseLong(tokens[2]);
				final long count_PASS = Long.parseLong(tokens[3]);
				Range r = this.ranges.stream().
						filter(R->R.lowerBound==m && R.upperBound==M).
						findFirst().
						orElse(null);
				if(r!=null) {
					r.count_ALL+= count_ALL;
					r.count_PASS+= count_PASS;
					}
				else {
					r=new Range(m,M);
					r.count_ALL = count_ALL;
					r.count_PASS = count_PASS;
					this.ranges.add(r);
					Collections.sort(this.ranges,(A,B)->Double.compare(A.lowerBound, B.lowerBound));
					}
				}
			}
		
		long getCountAll() {
			return ranges.stream().mapToLong(R->R.count_ALL).sum();
			}
		long getCountPass() {
			return ranges.stream().mapToLong(R->R.count_PASS).sum();
			}

		private OptionalDouble getPercentile(double f) {
			final long count = getCountAll();
			if(count==0L) return OptionalDouble.empty();
			long count2 = (long)(count*f);
			for(Range r: ranges) {
				for(long i=0;i< r.count_ALL;i++) {
					if(count2==0L) return OptionalDouble.of(r.lowerBound);
					count2--;
					}
				}
			throw new IllegalStateException();
			}
		
		OptionalDouble getLowPercentile() {
			return getPercentile(VcfGatkEval.this.percentile);
			}
		
		OptionalDouble getHighPercentile() {
			return getPercentile(1.0 - VcfGatkEval.this.percentile);
			}
		
		void saveAsR(PrintWriter out) throws IOException {
			if(ranges.isEmpty()) return;
			final long count_ALL = getCountAll();
			String title = getClass().getSimpleName();
			title = title.substring(0,title.length()-10);

			out.print("x <- c(");
			out.print(ranges.stream().map(R->String.valueOf((R.lowerBound+R.upperBound)/2.0)).collect(Collectors.joining(",")));
			out.println(")");
			
			out.print("y <- c(");
			out.print(ranges.stream().map(R->String.valueOf(R.count_ALL/(double)count_ALL)).collect(Collectors.joining(",")));
			out.println(")");
			
			double minx = 	this.ranges.stream().mapToDouble(R->R.lowerBound).min().getAsDouble();
			if(minx <=0 && logX()) minx=1.0;
			final double maxx = this.ranges.stream().mapToDouble(R->R.upperBound).max().getAsDouble() ;
			// plot ALL
			out.println("plot(x,y,pch=19,type=\"b\",col=\"blue\"," +
					"log=\""+(logX()?"x":"")+"\"," + 
					"xlim=c("+ minx +
						"," +  
						maxx +
						"),"+
					"ylim=c(0,max(y)),"+
					"sub=\"" + getSub() + "\"," +
					"main=\""+ title +" "+this.variantType.name()+"\",ylab=\"Density\",xlab=\""+getTag()+" "+this.variantType.name()+"\")");
			
			final long count_PASS = getCountPass();
			if(count_PASS>0L && 	this.ranges.stream().anyMatch(R->R.count_ALL!=R.count_PASS)) {
				out.print("y <- c(");
				out.print(ranges.stream().map(R->String.valueOf(R.count_PASS/(double)count_PASS)).collect(Collectors.joining(",")));
				out.println(")");
				// plot FILTERED
				out.println("lines(x,y,col=\"cyan4\",pch=19,type=\"b\")");

				out.println("legend(x=\"topright\",legend=c(\"ALL\",\"PASS\"),fill=c(\"blue\",\"cyan4\"))");

				}
			
			
			OptionalDouble limit = getLowPercentile();
			if(limit.isPresent()) {
				out.println("abline(v="+limit.getAsDouble()+", col=\"green\")");
				}
			limit = getHighPercentile();
			if(limit.isPresent()) {
				out.println("abline(v="+limit.getAsDouble()+", col=\"green\")");
				}
			}
		@Override
		public String toString() {
			return getTag()+"("+this.variantType.name()+")";
			}
	}
		
		private class QualByDepthAnnotation extends Annotation {
			QualByDepthAnnotation(VariantContext.Type t) { super(t);}
			@Override String getTag() { return GATKConstants.QD_KEY;}
			@Override int getPrecision() { return 1;}
			@Override OptionalDouble getHighPercentile() { return OptionalDouble.empty();}
			@Override String getSub() {return "Variant Confidence (QUAL) / Quality by Depth.";}

		}
		
		private class FisherStrandAnnotation extends Annotation {
			FisherStrandAnnotation(VariantContext.Type t) { super(t);}
			@Override String getTag() { return GATKConstants.FS_KEY;}
			@Override int getPrecision() { return 1;}
			@Override boolean logX() { return true;}
			@Override OptionalDouble getLowPercentile() { return OptionalDouble.empty();}
			@Override String getSub() {return "Phred-scaled p-value using Fisher's exact test to detect strand bias";}
		}
		
		private class StrandOddsRatioAnnotation extends Annotation {
			StrandOddsRatioAnnotation(VariantContext.Type t) { super(t);}
			@Override String getTag() { return GATKConstants.SOR_KEY;}
			@Override int getPrecision() { return 10;}
			@Override OptionalDouble getLowPercentile() { return OptionalDouble.empty();}
			@Override String getSub() {return "Symmetric Odds Ratio of 2x2 contingency table to detect strand bias";}
		}
		
		private class  RMSMappingQualityAnnotation extends Annotation {
			RMSMappingQualityAnnotation(VariantContext.Type t) { super(t);}
			@Override String getTag() { return GATKConstants.MQ_KEY;}
			@Override int getPrecision() { return 1;}
			@Override OptionalDouble getHighPercentile() { return OptionalDouble.empty();}
			@Override String getSub() {return "Mean square mapping quality over all the reads at the site";}
		}
		
		private class  MappingQualityRankSumTestAnnotation extends Annotation {
			MappingQualityRankSumTestAnnotation(VariantContext.Type t) { super(t);}
			@Override String getTag() { return GATKConstants.MQRankSum_KEY;}
			@Override int getPrecision() { return 10;}
			@Override OptionalDouble getHighPercentile() { return OptionalDouble.empty();}
			@Override String getSub() {return "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities";}

		}
		
		private class  DepthAnnotation extends Annotation {
			DepthAnnotation(VariantContext.Type t) { super(t);}
			@Override String getTag() { return VCFConstants.DEPTH_KEY;}
			@Override int getPrecision() { return 1;}
			@Override String getSub() {return "Depth";}
		}
		
		private class  ReadPosRankSumTestAnnotation extends Annotation {
			ReadPosRankSumTestAnnotation(VariantContext.Type t) { super(t);}
			@Override String getTag() { return "ReadPosRankSum";}
			@Override int getPrecision() { return 10;}
			@Override OptionalDouble getHighPercentile() { return OptionalDouble.empty();}
			@Override String getSub() {return "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias";}
		}



	public int doWork(List<String> args) {
		try {
			final List<Path> paths = IOUtils.unrollPaths(args);
			
		
			
			final VariantContext.Type[] ctx_types = new   VariantContext.Type[] {VariantContext.Type.SNP,VariantContext.Type.INDEL};
			final List<Annotation> annotations= new ArrayList<>();
			for(VariantContext.Type type:ctx_types) {
				annotations.add(new QualByDepthAnnotation(type));
				annotations.add(new FisherStrandAnnotation(type));
				annotations.add(new StrandOddsRatioAnnotation(type));
				annotations.add(new RMSMappingQualityAnnotation(type));
				annotations.add(new MappingQualityRankSumTestAnnotation(type));
				annotations.add(new ReadPosRankSumTestAnnotation(type));
				if(this.include_depth) {
					annotations.add(new DepthAnnotation(type));
					}
				}
			Collections.sort(annotations, (A,B)->{
				int i= A.getTag().compareTo(B.getTag());
				if(i!=0) return i;
				return A.variantType.compareTo(B.variantType);
				});
			
			switch (this.inputType) {
				case vcf: {
					if(paths.size()>1) {
						LOG.error("Expected zero or one vcf as input but got "+paths.size());
						return -1;
						}
					try(VCFIterator r = (paths.isEmpty()?
								new BcfIteratorBuilder().open(stdin()):
								new BcfIteratorBuilder().open(paths.get(0))
							)) {
						final VCFHeader header = r.getHeader();
						for(Annotation ann: annotations) {
							ann.init(header);
							}
						annotations.removeIf(A->A.infoHeaderLine==null);
						if(!annotations.isEmpty()) {
								while(r.hasNext()) {
									final VariantContext ctx = r.next();
									for(Annotation ann: annotations) {
										ann.visit(ctx);
									}
								}///end while
							}//end if
						}
					break;
					}
				case table : {
					final Consumer<BufferedReader> consummer = BR->{
						String line;
						try {
							while((line=BR.readLine())!=null) {
								final String tag = line;
								line = BR.readLine();
								if(line==null) throw new IOException("cannot read variant type");
								final VariantContext.Type type = VariantContext.Type.valueOf(line);
								final Annotation ann = annotations.stream().
										filter(A->A.getTag().equals(tag) && A.variantType.equals(type)).
										findFirst().orElse(null);
								if(ann==null) {
									throw new IOException("Cannot find annotation "+tag+":"+type);
									}
								ann.loadTable(BR);
								}
							}
						catch(IOException error ) {
							throw new RuntimeIOException(error);
							}
						};
					if(paths.isEmpty()) {
						try(BufferedReader br=super.openBufferedReader(null)) {
							consummer.accept(br);
							}
						}
					else
						{
						for(Path path:paths) {
							try(BufferedReader br= IOUtils.openPathForBufferedReading(path)) {
								consummer.accept(br);
								}
							}
						}
					break;
					}
				default: throw new IllegalStateException();
				}
			
			
		annotations.removeIf(A->A.getCountAll()==0L);
			
	
		try(PrintWriter w = IOUtils.openPathForPrintWriter(Paths.get(this.basename + ".output.R"))) {
			w.println("pdf(\"output.pdf\")");
			for(Annotation ann: annotations) {
					ann.saveAsR(w);
					}
			w.println("dev.off()");
			w.flush();
			}
		try(PrintWriter w = IOUtils.openPathForPrintWriter(Paths.get(this.basename + ".output.table.txt"))) {
			for(Annotation ann: annotations) {
					ann.saveTable(w);
					}
			w.flush();
			}
		
		try(PrintWriter w = IOUtils.openPathForPrintWriter(Paths.get(this.basename + ".output.filters.txt"))) {
			for(Annotation ann: annotations) {
					ann.saveFilter(w);
					}
			w.flush();
			}
		return 0;
	    }
	catch(final Throwable err ) {
		LOG.error(err);
		return -1;
	    }
    }

	public static void main(String[] args) {
		new VcfGatkEval().instanceMainWithExit(args);
	}

}
