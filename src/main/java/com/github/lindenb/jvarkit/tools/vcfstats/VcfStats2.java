/*
The MIT License (MIT)
Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.DoubleConsumer;
import java.util.function.DoubleUnaryOperator;

import org.springframework.context.ApplicationContext;
import org.springframework.context.support.ClassPathXmlApplicationContext;
import org.springframework.context.support.FileSystemXmlApplicationContext;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## Example

```
java -jar dist/jvarkit.jar vcfstats src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz |  R --no-save 
```

END_DOC
 */
@Program(name="vcfstats2",
	description="Produce VCF statitics",
	keywords={"vcf","stats","multiqc"},
	modificationDate = "20240726",
	creationDate = "20131212",
	jvarkit_amalgamion = true,
	jvarkit_hidden = true,
	generate_doc = false,
	menu="VCF Manipulation"
	)
public class VcfStats2 extends Launcher {
	private static final Logger LOG = Logger.of(VcfStats2.class);
	public static final String DEFAULT_MAIN_BEAN_NAME="main";


	@Parameter(names={"-o","--output"},description="output directory",required=true)
	private Path outputDirectory = null;
	
	@Parameter(names={"--pipe"},description="write input VCF to stdout")
	private boolean enabled_pipe_to_stdout = false;

	@Parameter(names={"--prefix"},description="prefix for output files")
	private String prefix = "vcfstats.";

	
	@Parameter(names={"--title"},description="main section title")
	private String main_section_title = "";
	@Parameter(names={"--description"},description="main section description")
	private String main_section_description = "";
	
	@Parameter(names={"--sample2population","--sample2pop"},description=SamplePopulation.OPT_DESC)
	private Path sample2catPath = null;
	@Parameter(names={"--other-samples"},description="if sample doesn't belong to a population in sample2pop, "
			+ "create a new population named 'x' and insert those lonely samples in that population")
	private String otherPopulationName="other";

	
	@Parameter(names={"-exclude","--exclude"},description="name of modules to be excluded",hidden=true)
	private String moduleExcludeStr = "";
	@Parameter(names={"--list"},description="list available modules and exit",help = true)
	private boolean list_modules = false;	
	
	@Parameter(names={"--xml","--spring"},description="XML Spring config")
	private String xmlSpringConfigLocation=null;
	@Parameter(names={"--bean"},description="main bean id in the XML configuration")
	private String mainBeanId = DEFAULT_MAIN_BEAN_NAME;
	
	private final SamplePopulation samplePopulation=new SamplePopulation();
	
	/** manage x =  ((int)(x/factor))*factor   */
	private class DoubleRounder implements DoubleUnaryOperator {
		final DecimalFormat decimalFormat;
		DoubleRounder(int n) {
			String str="#.";
			while(n>0) {
				str+="#";
				--n;
				}
			this.decimalFormat = new DecimalFormat(str);
			this.decimalFormat.setRoundingMode(java.math.RoundingMode.CEILING);
			}
		@Override
		public double applyAsDouble(double v) {
			return Double.parseDouble(this.decimalFormat.format(v));
			}
		}
	
	private class BinRounder implements DoubleUnaryOperator {
		final double[] bins;
		BinRounder(double...values) {
			bins = Arrays.copyOf(values, values.length);
			Arrays.sort(bins);
			if(bins.length==0) throw new IllegalArgumentException();
			}
		@Override
		public double applyAsDouble(double v) {
			for(int i=0;i< this.bins.length;i++) {
				if(v <= this.bins[i]) return this.bins[i];
				}
			return this.bins[this.bins.length-1];
			}
		}
	
	private static class DataPoint implements DoubleConsumer {
		private long count=0L;
		private double sum = 0.0;
		
		@Override
		public void accept(double value) {
			count++;
			sum+=value;
			}
		
		public long getCount() {
			return this.count;
			}
	
		public  OptionalDouble getAverage() {
			return isEmpty()?OptionalDouble.empty():OptionalDouble.of(this.sum/this.count);
			}
		public boolean isEmpty() {
			return this.getCount()==0;
			}
		}

	private void loadPhenotypes(final VCFHeader header) throws IOException {
		if(!header.hasGenotypingData()) return;
		final Set<String> sampleset = new HashSet<>(header.getGenotypeSamples());
		if(sample2catPath!=null) {
			samplePopulation.load(sample2catPath);
			samplePopulation.retain(header);
			}
		if(!StringUtils.isBlank(otherPopulationName)) {
			if(this.samplePopulation.hasCollection(this.otherPopulationName)) {
				throw new IllegalArgumentException("collection \""+this.otherPopulationName+"\" already defined in "+this.sample2catPath);
				}
			for(String sn:sampleset) {
				if(samplePopulation.hasSample(sn)) continue;
				samplePopulation.insert(sn, this.otherPopulationName);
				}
			}
		}
	
	
	
	@Override
	public int doWork(final List<String> args) {
			
		try {
			final ApplicationContext springApplicationContext;
			if(StringUtils.isBlank(this.xmlSpringConfigLocation))
				{
				this.xmlSpringConfigLocation = "META-INF/spring/default-vcfstats-analyzers.xml";
				LOG.warn("No spring xml was defined. Using default in jar file jar://"+this.xmlSpringConfigLocation);
				springApplicationContext = new ClassPathXmlApplicationContext(this.xmlSpringConfigLocation);
				}
			else
				{
				springApplicationContext =  new FileSystemXmlApplicationContext(xmlSpringConfigLocation);
				}
		
			if(!springApplicationContext.containsBean(this.mainBeanId))
				{
				LOG.error("cannot get bean "+this.mainBeanId+" in "+ this.xmlSpringConfigLocation);
				return -1;
				}
			
			final Object o = springApplicationContext.getBean(this.mainBeanId);
			if( o == null) {
				throw new IllegalArgumentException("bean "+ this.mainBeanId +" was found but it is null in " + this.xmlSpringConfigLocation);
				}
			if( !(o instanceof List)) {
				throw new IllegalArgumentException("bean "+ this.mainBeanId +" is not an instance of java.util.List in " + this.xmlSpringConfigLocation);
				}
			
			final List<?> list = (List<?>)o;
			final List<Analyzer> modules =new ArrayList<>(list.size());
			
			for(int i=0;i< list.size();i++) {
				final Object o2 = list.get(i);
				if( o2 == null) {
					throw new IllegalArgumentException("element["+i+"]found but it is null in "+ this.xmlSpringConfigLocation);
					}
				if( !(o2 instanceof Analyzer)) {
					throw new IllegalArgumentException("element["+i+"]found but it not an instance of "+Analyzer.class.getName()+" in "+   this.xmlSpringConfigLocation);
					}
				final Analyzer analyzer = Analyzer.class.cast(o2);
				if(modules.stream().anyMatch(P->P.getName().equals(analyzer.getName()))) {
					LOG.error("duplicate analyzer with name "+analyzer.getName());
					return -1;
					}
				modules.add(Analyzer.class.cast(analyzer));
				}
			
	
	
	// remove modules
	modules.removeIf(M->Arrays.stream(this.moduleExcludeStr.split("[,; \t:]")).anyMatch(S->S.equalsIgnoreCase(M.getName())));
		

	if(this.list_modules) {
		try(PrintWriter w = super.openPathOrStdoutAsPrintWriter(null)) {
			for(Analyzer analyzer:modules) {
				w.print(analyzer.getName());
				w.print("\t");
				w.println(analyzer.getDescription());
				}	
			w.flush();
			}
		return 0;
		}
	
	final String input = oneFileOrNull(args);
	IOUtil.assertDirectoryIsWritable(this.outputDirectory);
	
	final VariantContextWriter variantCtxWriter =  enabled_pipe_to_stdout?
			VCFUtils.createVariantContextWriterToStdout():
			null;
	
	try(VCFIterator iter= super.openVCFIterator(input)) {
		final VCFHeader header=iter.getHeader();
		loadPhenotypes(header);
		if(variantCtxWriter!=null) {
			variantCtxWriter.writeHeader(header);
			}
		
		for(Analyzer analyzer:modules) {
			analyzer.init(header,this.samplePopulation);
			if(!analyzer.isEnabled()) {
				LOG.warn("module "+analyzer.getName()+" will be disabled. ["+analyzer.getClass().getSimpleName()+"]");
				}
			}
		modules.removeIf(M->!M.isEnabled());
		if(modules.isEmpty()) {
			LOG.warn("no module was enabled");
			}
		while(iter.hasNext()) {
			final VariantContext ctx = iter.next();
			if(variantCtxWriter!=null) variantCtxWriter.add(ctx);
			for(Analyzer analyzer:modules) {
				analyzer.visit(ctx);
				}
			}
		if(variantCtxWriter!=null) {
			variantCtxWriter.close();
			}
		
		for(Analyzer analyzer:modules) {
			analyzer.writeReports(this.outputDirectory,FunctionalMap.make());
			}
		}
		
		return 0;
		}
	catch (final Throwable e) {
		e.printStackTrace();
		LOG.error(e);
		return -1;
		}

	}
	
			
	public static void main(final String[] args)
		{
		new VcfStats2().instanceMainWithExit(args);
		}
	}
