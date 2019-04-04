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
package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.Vector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.chart.BarChart;
import com.github.lindenb.jvarkit.chart.CategoryAxis;
import com.github.lindenb.jvarkit.chart.Chart;
import com.github.lindenb.jvarkit.chart.LineChart;
import com.github.lindenb.jvarkit.chart.NumberAxis;
import com.github.lindenb.jvarkit.chart.RExporter;
import com.github.lindenb.jvarkit.chart.StackedBarChart;
import com.github.lindenb.jvarkit.chart.XYChart;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.RangeOfIntegers;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/*
BEGIN_DOC

## Examples

todo

## Screenshot

todo

END_DOC
*/
@Program(name="bamstatsjfx",
description="BAM statistics",
keywords={"bam","stats"},
generate_doc=false
)
public class BamStatsJfx extends Launcher {
	private static final Logger LOG=Logger.build(BamStatsJfx.class).make();

	@Parameter(names={"-s","--seconds"},description="Save every 's' seconds if output file is defined")
	private int refreshEverySeconds = 15;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--stdin"},description="if there is no file argument. Read vcf from stdin instead of opening a FileOpen dialog")
	private boolean vcf_stdin = false;
	@Parameter(names={"-hr","--hr"},description="Ignore HOM_REF in genotype type.")
	private boolean ignore_HOM_REF = false;
	@Parameter(names={"-fgt","--fgt"},description="Ignore filtered **GENOTYPES**")
	private boolean ignore_filtered_genotypes = false;
	@Parameter(names={"--prefix"},description="Title Prefix")
	private String titlePrefix="";
	@Parameter(names={"--read-lengths"},description="tranches for Read lengths. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers readLengthTranches = new RangeOfIntegers(0,10,20,30,40,50,100,150,200,300,400,500,1000);
	@Parameter(names={"--max-read-length"},description="max read length for base composition. 0: ignore; <0 not limit")
	private int max_read_length = 500;
	@Parameter(names={"--clip-lengths"},description="tranches for Clip lengths. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers clipLengthTranches = new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50);

	
	
	
	private abstract class ChartGenerator
		{
		boolean enabled = true;
		long nRecords = 0;
		boolean isEnabled() {
			return enabled;
			}
		String getChartTitle() {
			return "untitled";
		}
			
		abstract Chart makeChart();
		
		void title(final Chart c,final String s) {
			c.setTitle(
					StringUtil.isBlank(titlePrefix)?
					s:
					titlePrefix+ " : " + s
					);
			}
		
		abstract void visit(final SAMRecord ctx) ;
		
		// public String getFilename() { return getChartTitle().replaceAll("[^A-Za-z_0-9]+","")+".png";}
		}
	
	private class ReadLengthGenerator extends ChartGenerator
		{
		private final Counter<RangeOfIntegers.Range> len2count = new Counter<>();
		@Override
		String getChartTitle() {
			return "Read-Length";
			}
		
		
		@Override
		Chart makeChart() {

			if(this.len2count.isEmpty()) return null;
			
			final List<String> lengthCategories = this.len2count.keySet().
					stream().
					sorted().
					map(x->x.toString()).
					collect(Collectors.toList());
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis( lengthCategories);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			
			final BarChart<String,Number> bc = 
		            new BarChart<>(xAxis,yAxis);
		
			
			for(final RangeOfIntegers.Range x:new TreeSet<>(this.len2count.keySet()))
				{
				
				series1.getData().add(new XYChart.Data<String,Number>(
						x.toString(),
						this.len2count.count(x)
					));
						
				}
			bc.getData().add(series1);
			bc.setCategoryGap(1);
			title(bc,
	        		this.getChartTitle()+
	        		" N. Primary SAMRecord "+ StringUtils.niceInt(nRecords)
	        		)
	        		;
	        xAxis.setLabel("Length");
	        yAxis.setLabel("Number of Reads");
	        xAxis.setTickLabelRotation(90);
	        bc.setLegendVisible(false);
	        return bc;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			this.nRecords ++ ;
			final RangeOfIntegers.Range r = readLengthTranches.getRange(rec.getReadLength());
			this.len2count.incr(r);
			}
		
		}
	private class ClippedReadsGenerator extends ChartGenerator
		{
		private final Counter<RangeOfIntegers.Range> clip5 = new Counter<>();
		private final Counter<RangeOfIntegers.Range> clip3 = new Counter<>();
		
		@Override
		String getChartTitle() {
			return "Clip";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(clipLengthTranches.getRanges().stream().map(O->O.toString()).collect(Collectors.toList()));
			
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
			
	
			for(int side=0;side<2;++side)
				{
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(side==0?"5'":"3'");
				for(final RangeOfIntegers.Range r: clipLengthTranches.getRanges()) {
					series1.getData().add(new XYChart.Data<String,Number>(
							r.toString(),
							(side==0?this.clip5:this.clip3).count(r)
							));
					}
				bc.getData().add(series1);
				}
			bc.setCategoryGap(1);
			bc.setLegendVisible(true);;
			title(bc,this.getChartTitle()+" N. Clipped Records: "+StringUtils.niceInt(nRecords));
	        xAxis.setLabel("Clip-Length");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			if(rec.getReadUnmappedFlag()) return ;
			final Cigar cigar = rec.getCigar();
			if(cigar==null || cigar.isEmpty()) return;
			if(!cigar.isClipped()) return;
			this.nRecords++;
			CigarElement ce =cigar.getFirstCigarElement();
			if(ce.getOperator().isClipping()) {
				this.clip5.incr(clipLengthTranches.getRange(ce.getLength()));
			}
			
			ce =cigar.getLastCigarElement();
			if(ce.getOperator().isClipping()) {
					this.clip3.incr(clipLengthTranches.getRange(ce.getLength()));
				}
			}
		}
	
	private class BaseCompositionGenerator extends ChartGenerator
		{
		private final List<Counter<Character>> pos2count = new ArrayList<>();
		BaseCompositionGenerator() {
			}
		
		@Override
		String getChartTitle() {
			return "Base Composition";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(
						IntStream.
						range(0, pos2count.size()).
						mapToObj(x->StringUtils.niceInt(x+1)).
						collect(Collectors.toList())
						);
			final TreeSet<Character> bases = pos2count.stream().
					flatMap(F->F.keySet().stream()).
					collect(Collectors.toCollection(TreeSet::new));
			
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
		

			for(final Character base:bases)
				{
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(base.toString());
				for(int x=0;x< this.pos2count.size();++x) {
					series1.getData().add(new XYChart.Data<String,Number>(
							StringUtils.niceInt(x+1),
							this.pos2count.get(x).count(base)
							));
					}
				bc.getData().add(series1);
				}
			bc.setCategoryGap(1);
			title(bc,this.getChartTitle()+" N. Records: "+StringUtils.niceInt(nRecords));
	        xAxis.setLabel("Position");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Bases");
	        return bc;
			}
		
		private void visit(char c,int pos) {
			this.pos2count.get(pos).incr(Character.toUpperCase(c));
		}
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			this.nRecords++;
			final byte[] readBases = rec.getReadBases();
			
			while(this.pos2count.size()<=readBases.length && (max_read_length<0 || this.pos2count.size()<max_read_length)) {
				this.pos2count.add(new Counter<>());
				}
			
			if(rec.getReadNegativeStrandFlag())
				{
				for(int i=0;i< readBases.length && (max_read_length<0 || i<max_read_length);i++)
					{
					visit(AcidNucleics.complement((char)readBases[(readBases.length-1)-i]),i);
					}
				}
			else
				{
				for(int i=0;i< readBases.length && (max_read_length<0 || i<max_read_length);i++)
					{
					visit((char)readBases[i],i);
					}
				}
			}
		}

	private class QualGenerator extends ChartGenerator
		{
		private class Average {
			double sum=0.0;
			long n=0L;
			double get() {
				return n==0?0.0:sum/n;
			}
		}
		private final Map<Integer,Average> pos2qual = new HashMap<>();
		
		@Override
		String getChartTitle() {
			return "Quality";
			}
		
		@Override
		Chart makeChart() {

			if(this.pos2qual.isEmpty()) return null;
			
			final int max_length = pos2qual.keySet().stream().mapToInt(X->X.intValue()).max().getAsInt();
			
			final NumberAxis yAxis = new NumberAxis();
			final NumberAxis xAxis = new NumberAxis();
			
			final LineChart<Number,Number> bc = 
		            new LineChart<Number,Number>(xAxis,yAxis);
			final XYChart.Series<Number,Number> series1 = new XYChart.Series<>();
			for(int x=0;x<= max_length;++x)
				{
				series1.getData().add(new XYChart.Data<Number,Number>(
					x+1,
					pos2qual.get(x).get()
					));
				}
			bc.getData().add(series1);
						
			bc.setLegendVisible(false);
			title(bc,this.getChartTitle());
	        xAxis.setLabel("Position");
			xAxis.setLowerBound(0.0);
			xAxis.setUpperBound((double)max_length);
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Avg Qual.");
	        yAxis.setLowerBound(0.0);
	        return bc;
			}
		
		private void visit(int x,byte q) {
			Average a = this.pos2qual.get(x);
			if(a==null) {
				a = new Average();
				this.pos2qual.put(x, a);
				}
			a.n++;
			a.sum+=q;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			final byte[] readQuals = rec.getBaseQualities();
			if(readQuals.length==0) return;
			this.nRecords++;
			for(int i=0;i< readQuals.length ;i++)
				{
				int j = rec.getReadNegativeStrandFlag()? (readQuals.length-1-i):i;
				visit(j,readQuals[j]);	
				}
			}
		}
	
	private class GCPercentGenerator extends ChartGenerator
		{
		private final Counter<Integer> mappedgc2count = new Counter<>();
		private final Counter<Integer> unmappedgc2count = new Counter<>();
		GCPercentGenerator() {
			}
		
		@Override
		String getChartTitle() {
			return "GC Percent";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final NumberAxis xAxis = new NumberAxis();
			
			final LineChart<Number,Number> bc = 
		            new LineChart<Number,Number>(xAxis,yAxis);
	        for(int side=0;side< 2;++side) {
				final XYChart.Series<Number,Number> series1 = new XYChart.Series<>();
				series1.setName(side==0?"Mapped":"Unmapped");
				final Counter<Integer> rc =(side==0?this.mappedgc2count:this.unmappedgc2count);
				for(final Integer percent: new TreeSet<>(rc.keySet()))
					{
					series1.getData().add(new XYChart.Data<Number,Number>(
						percent,
						rc.count(percent)
						));
					}
				bc.getData().add(series1);
				}			
			bc.setLegendVisible(true);
			title(bc,this.getChartTitle());
	        xAxis.setLabel("GC%");
			xAxis.setLowerBound(0.0);
			xAxis.setUpperBound(100.0);
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Count");
	        yAxis.setLowerBound(0.0);
	        return bc;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			this.nRecords++;
			final byte[] readBases = rec.getReadBases();
			if(readBases.length==0) return;
			double ngc=0;
			for(int i=0;i< readBases.length;++i) {
				switch(readBases[i])
					{
					case 's':case 'S':
					case 'g':case 'G':
					case 'c':case 'C': ngc++;break;
					}
				}
			final int ngi = (int)((ngc/(double)readBases.length)*100.0);
			if(rec.getReadUnmappedFlag())
				{
				this.unmappedgc2count.incr(ngi);
				}
			else
				{
				this.mappedgc2count.incr(ngi);
				}
			}
		}
	
	
	private class ContigUsageGenerator extends ChartGenerator
		{
		private final Counter<String> count = new Counter<>();
		private final SAMSequenceDictionary dict;
		ContigUsageGenerator(final SAMSequenceDictionary dict) {
			this.dict = dict;
			}
		
		@Override
		String getChartTitle() {
			return "Contigs";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final List<String> contigs = 
					this.dict.getSequences().stream().
						map(S->S.getSequenceName()).
						filter(S->count.count(S)>0).
						collect(Collectors.toList())
						;
			if(this.count.count(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)>0L)
				{
				contigs.add(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
				}
			
			final CategoryAxis xAxis = new CategoryAxis(contigs);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final SAMSequenceRecord ssr: this.dict.getSequences())
				{
				long n = this.count.count(ssr.getSequenceName());
				if(n==0L) continue;
				series1.getData().add(new XYChart.Data<String,Number>(ssr.getSequenceName(),n));
				}
			long n=this.count.count(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
			if(n>0L)
				{
				series1.getData().add(new XYChart.Data<String,Number>(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME,n));
				}
			bc.getData().add(series1);
			title(bc,this.getChartTitle()+ " Num. Reads. "+StringUtils.niceInt(this.nRecords));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Contig");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			this.nRecords++;
			if(rec.getReadUnmappedFlag())
				{
				this.count.incr(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
				}
			else
				{
				this.count.incr(rec.getContig());
				}
			}
		}
	
	private class SamFlagUsagGenerator extends ChartGenerator
		{
		private final Counter<SAMFlag> count = new Counter<>();
		
		@Override
		String getChartTitle() {
			return "SamFlags";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final List<String> contigs = Arrays.asList(SAMFlag.values()).stream().map(S->S.name()).collect(Collectors.toList());
			
			final CategoryAxis xAxis = new CategoryAxis(contigs);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final SAMFlag f: SAMFlag.values())
				{
				series1.getData().add(new XYChart.Data<String,Number>(f.name(),count.count(f)));
				}
			bc.getData().add(series1);
			title(bc,this.getChartTitle()+ " Num. Reads. "+ StringUtils.niceInt(this.nRecords));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Flag");
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			this.nRecords++;
			for(final SAMFlag f: rec.getSAMFlags())
				{
				count.incr(f);
				}
			}
		}
	
	private final List<ChartGenerator> chartGenerators = new Vector<>();
	
	
	 private void save() throws IOException {
			try (final PrintWriter pw = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
			final RExporter exporter = new RExporter();
			for(final ChartGenerator cg: this.chartGenerators) {
				LOG.info("saving "+cg.getChartTitle());
				if(!cg.isEnabled()) continue;
				final Chart chart = cg.makeChart();
				if(chart==null) continue;
				exporter.exportToR(pw, chart);
				}
			pw.flush();
			}
		}
	 
	@Override
	public int doWork(final List<String> args) {
		SamReader samReader = null;
		SAMRecordIterator iter = null;
		try {
			samReader = super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader header =  samReader.getFileHeader();
			
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict!=null && !dict.isEmpty()) {
				chartGenerators.add(new ContigUsageGenerator(dict));
				}
			this.chartGenerators.add(new ReadLengthGenerator());
			if(max_read_length!=0) { /* negative : not limit */
				this.chartGenerators.add(new BaseCompositionGenerator());
				}
			this.chartGenerators.add(new GCPercentGenerator());
			this.chartGenerators.add(new QualGenerator());
			this.chartGenerators.add(new SamFlagUsagGenerator());
			this.chartGenerators.add(new ClippedReadsGenerator());
			
			
			
			
			
			
			chartGenerators.removeIf(G->!G.isEnabled());
			
			iter = samReader.iterator();
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			long last = -1L;
			while(iter.hasNext())
				{
				final SAMRecord rec = progress.apply(iter.next());
				
				for(final ChartGenerator cg: this.chartGenerators) {
					if(!cg.isEnabled()) continue;
					cg.visit(rec);
					}
				if(this.outputFile!=null)
					{
					final long now = System.currentTimeMillis();
					if(last<0L || now - last > this.refreshEverySeconds * 1000) {
						last=now;
						save();
						}
					}
				}
			progress.close();
			iter.close();iter=null;
			samReader.close();samReader=null;
			save();
			return 0;
			}
		catch(final Exception err) {
			err.printStackTrace();
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			}
		}
	
	

	public static void main(final String[] args) {
		new BamStatsJfx().instanceMainWithExit(args);
	}

}
