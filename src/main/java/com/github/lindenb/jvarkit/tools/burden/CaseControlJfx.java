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
package com.github.lindenb.jvarkit.tools.burden;


import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.chart.NumberAxis;
import com.github.lindenb.jvarkit.chart.RExporter;
import com.github.lindenb.jvarkit.chart.ScatterChart;
import com.github.lindenb.jvarkit.chart.XYChart;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory.AFExtractor;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;


/**
BEGIN_DOC

### Example

```
java -jar dist/casectrljfx.jar --pedigree  mutations.ped mutations.vcf
```
### See also:

  * https://twitter.com/yokofakun/status/860495863633805312

## screenshot

![screenshot](https://pbs.twimg.com/media/C_EYa54W0AAopkl.jpg)

## History

  * removed JFX because openjdk doesn't support JFX

END_DOC
 */
@Program(
	name="casectrljfx",
	description="chart of case/control maf from a VCF and a pedigree",
	keywords={"vcf","pedigree","case","control","visualization","jfx","chart","maf"}
	)
public class CaseControlJfx extends Launcher {
	
	private static final Logger LOG = Logger.build(CaseControlJfx.class).make();

	private enum VariantPartitionType{
		chromosome,
		variantType,
		autosomes,
		qual,
		vqslod,
		typeFilter,
		distance,
		n_alts
		}
	
	private static interface VariantPartition
		{
		public List<XYChart.Series<Number,Number>> getSeries();
		public void add(VariantContext vc,Pedigree ped, XYChart.Data<Number,Number> data);
		}
	
	private static class VariantTypePartition implements VariantPartition
		{
		final List<XYChart.Series<Number,Number>> series = new ArrayList<>( VariantContext.Type.values().length);
		VariantTypePartition() 
			{
			for(final VariantContext.Type t:VariantContext.Type.values()) {
				final XYChart.Series<Number,Number> serie = new XYChart.Series<>();
				serie.setName(t.name());
				this.series.add(serie);
				}
			}
		@Override
		public List<XYChart.Series<Number,Number>> getSeries()
			{
			return this.series;
			}
		@Override
		public void add(VariantContext vc,Pedigree ped,XYChart.Data<Number,Number> data){
			if(vc==null) return;
			final VariantContext.Type t= vc.getType();
			this.series.get(t.ordinal()).getData().add(data);
			}
		}
	
	private static abstract class MapPartition implements VariantPartition
		{
		final Map<String,XYChart.Series<Number,Number>> series = new TreeMap<>();
		MapPartition() 
			{
			}
		
		protected XYChart.Series<Number,Number> get(String key) {
			XYChart.Series<Number,Number> S = this.series.get(key);
			if(S==null)
				{
				S = new XYChart.Series<>();
				S.setName(key);
				this.series.put(key,S);
				}
			return S;
			}
		@Override
		public List<XYChart.Series<Number,Number>> getSeries()
			{
			return new ArrayList<>(this.series.values());
			}
	
		}


	private static class ChromosomePartition extends MapPartition
		{
		ChromosomePartition() 
			{
			}
		
		protected String normalizeContig(final String str) {
			return str;
			}
		
		@Override
		public void add(VariantContext vc,Pedigree ped,XYChart.Data<Number,Number> data){
			if(vc==null) return;
			final String str = normalizeContig(vc.getContig());
			if(str==null || str.isEmpty()) return;
			this.get(str).getData().add(data);
			}
		}
	
	private static class SexualContigPartition extends ChromosomePartition
		{
		protected String normalizeContig(final String str) {
			if( str.equalsIgnoreCase("chrX") || str.equalsIgnoreCase("X") || 
				str.equalsIgnoreCase("chrY") || str.equalsIgnoreCase("Y")
				) return str; 
			return "Autosome";
			}
		}

	private static class QualPartition extends MapPartition
		{
		@Override
		public void add(VariantContext vc,Pedigree ped,XYChart.Data<Number,Number> data){
			if(vc==null) return;
			final String str ;
			if(vc.hasLog10PError())
				{		
				int qual=(int)vc.getPhredScaledQual();
				if(qual <10) str="<10";
				else if(qual <20) str="<20";
				else if(qual <30) str="<30";
				else if(qual <40) str="<40";
				else if(qual <50) str="<50";
				else if(qual <100) str="<100";
				else str=">=100";
				}
			else
				{
				str="N/A";
				}
			this.get(str).getData().add(data);
			}
		}
	private static class VQSLODPartition extends MapPartition
		{
		@Override
		public void add(VariantContext vc,Pedigree ped,XYChart.Data<Number,Number> data){
			if(vc==null) return;
			if(!vc.hasAttribute("VQSLOD")) return;
			double vqslod=vc.getAttributeAsDouble("VQSLOD",Double.NaN);
			if(Double.isNaN(vqslod)) return;
			final int WINDOW=5;
			int norm=WINDOW*(int)(vqslod/((double)WINDOW));
			String str=String.valueOf(norm);
			
			this.get(str).getData().add(data);
			}
		}

	private static class NAltsPartition extends MapPartition
		{
		@Override
		public void add(VariantContext vc,Pedigree ped,XYChart.Data<Number,Number> data){
			if(vc==null) return;
			this.get(String.valueOf(vc.getAlternateAlleles().size())).getData().add(data);
			}
		}
	
	private static class DisanceToDiagonalPartiton extends MapPartition
		{
		@Override
		public void add(VariantContext vc,Pedigree ped,XYChart.Data<Number,Number> data){
			if(vc==null) return;
			
			double x1= data.getXValue().doubleValue();
			double y1= data.getYValue().doubleValue();
			// use 'mirror' point on diagonal
			double x2 = y1;
			double y2 = x1;
			//distance to diagonal is half distance between (x1,y1) and (x2,y2)
			double distance= Math.sqrt( Math.pow((x1-x2),2)+ Math.pow(y1-y2,2))/2.0;
			final double WINDOW=10.0;
			double norm= Math.abs(((int)(distance*WINDOW))/(double)WINDOW);
			this.get(String.valueOf(norm)).getData().add(data);
			}
		}
	
	private static class TypeAndFilterPartiton extends MapPartition
		{
		@Override
		public void add(VariantContext vc,Pedigree ped,XYChart.Data<Number,Number> data){
			if(vc==null) return;
			
			final String str=vc.getType().name()+(vc.isFiltered()?" FILTERED":" PASS");
			
			this.get(str).getData().add(data);
			}
		}
	
	
	private enum SelectSamples { all,males,females};
	
	
		@Parameter(names={"-p","--ped","--pedigree"},description="Pedigree File. If not defined, try to use the pedigree inserted in the VCF header.")
		File pedigreeFile = null;
		@Parameter(names={"-partition","--partition"},description="partition type. How series are built. For example 'variantType' will produces some series for INDEL, SNP, etc... ")
		VariantPartitionType partitionType=VariantPartitionType.variantType;
		@Parameter(names={"-nchr","--nocallishomref"},description="treat no call as HomRef")
		boolean no_call_is_homref=false;
		@Parameter(names={"-filter","--filter"},description="Ignore FILTERed variants")
		boolean ignore_ctx_filtered=false;
		@Parameter(names={"-gtfilter","--genotypefilter"},description="Ignore FILTERed Genotypes")
		boolean ignore_gt_filtered=false;
		//@Parameter(names={"--legendside"},description="Legend side")
		//Side legendSide = Side.RIGHT;
		//@Parameter(names={"--tooltip"},description="add mouse Tooltip the point (requires more memory)")
		//boolean add_tooltip = false;
		@Parameter(names={"--limit"},description="Limit to 'N' variants. negative==no limit; All point are loaded in memory. The more variants you have, the more your need memory")
		int limit_to_N_variants = -1;
		@Parameter(names={"--sex"},description="Select/Filter samples on their gender.")
		SelectSamples selectSamples=SelectSamples.all;
		@Parameter(names={"--title"},description="Default title for the graph")
		String userTitle=null;
		@Parameter(names={"--opacity"},description="Point opaciy [0-1]")
		double dataOpacity=0.4;
		@Parameter(names={"-o","--out"},description="Save the image in a file and then exit.")
		File outputFile=null;
		@Parameter(names={"-mafTag","--mafTag"},description=
				"[20180905] Do not calculate MAF for controls, but use this tag to get Controls' MAF. " +
				AFExtractorFactory.OPT_DESC
				)
		String controlFields =null;
		
		public CaseControlJfx()
			{
			}
		
		@Override
		public int doWork(final List<String> args) {
			final VariantPartition partition;
			Pedigree pedigree = null;
			VCFIterator in = null;
			try {
				
				switch(this.partitionType)
					{
					case variantType: partition = new VariantTypePartition(); break;
					case chromosome: partition = new ChromosomePartition(); break;
					case autosomes: partition = new SexualContigPartition(); break;
					case qual : partition = new QualPartition();break;
					case vqslod : partition = new VQSLODPartition();break;
					case typeFilter : partition = new TypeAndFilterPartiton();break;
					case distance : partition = new DisanceToDiagonalPartiton(); break;
					case n_alts : partition = new NAltsPartition(); break;
					default: throw new IllegalStateException(this.partitionType.name());
					}
				in = openVCFIterator(oneFileOrNull(args));
				
				if(this.pedigreeFile!=null)
					{
					pedigree = Pedigree.newParser().parse(this.pedigreeFile);
					}
				else
					{
					pedigree = Pedigree.newParser().parse(in.getHeader());
					}
				final AFExtractor controlAFExtractor;
				if(!StringUtil.isBlank(this.controlFields))
					{
					final List<AFExtractor> extractors = new AFExtractorFactory().parseFieldExtractors(this.controlFields);
					if(extractors.size()!=1) {
						LOG.error("extractor list should have size==1 . got "+extractors);
						return -1;
						}
					controlAFExtractor = extractors.get(0);
					if(!controlAFExtractor.validateHeader(in.getHeader())) {
						LOG.error("Invalid : "+controlAFExtractor);
						return -1;
						}
					}
				else
					{
					controlAFExtractor = null;
					}
				
				int count = 0;
				final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(in.getHeader()).logger(LOG).build();
				while(in.hasNext() && (this.limit_to_N_variants<0 || count<this.limit_to_N_variants)) 
					{
					final VariantContext ctx=progress.apply(in.next());
					
					if(this.ignore_ctx_filtered && ctx.isFiltered()) continue;
					
					++count;
					
					final List<Allele> alternates = ctx.getAlternateAlleles();
					for(int alt_idx=0;alt_idx < alternates.size();++alt_idx) {
						final Allele alt = alternates.get(alt_idx);
						final Double mafs[]={null,null};
						
						for(int i=0;i< 2;++i)
							{
							if(i==1 && controlAFExtractor!=null)
								{
								final List<Double> dvals = controlAFExtractor.parse(ctx);
								if(alt_idx< dvals.size() && dvals.get(alt_idx)!=null) {	
								    final double d= dvals.get(alt_idx);
									if(!Double.isNaN(d) && d>=0 && d<=1.0) mafs[1]=d;
									}
								}
							else
								{
								final MafCalculator mafCalculator = new MafCalculator(alt, ctx.getContig());
								mafCalculator.setNoCallIsHomRef(no_call_is_homref);
								for(Pedigree.Person person: (i==0?pedigree.getAffected():pedigree.getUnaffected()))
									{
									if(selectSamples.equals(SelectSamples.males) && !person.isMale()) continue;
									if(selectSamples.equals(SelectSamples.females) && !person.isFemale()) continue;
									
									final Genotype genotype = ctx.getGenotype(person.getId());
									if(genotype==null) continue;
									if(ignore_gt_filtered && genotype.isFiltered()) continue;
									mafCalculator.add(genotype, person.isMale());
									}
								if(!mafCalculator.isEmpty())
									{
									mafs[i]=mafCalculator.getMaf();
									}
								}
							}
						if(mafs[0]==null || mafs[1]==null) continue;
						final XYChart.Data<Number,Number> data = new  XYChart.Data<Number,Number>(mafs[0],mafs[1]);
						partition.add(ctx,pedigree,data);
						}
					}
				progress.close();
				in.close();in=null;
				
				
		        final NumberAxis xAxis = new NumberAxis(0.0,1.0,0.1);
		        xAxis.setLabel("Cases");
		        
		        final NumberAxis yAxis = new NumberAxis(0.0,1.0,0.1);
		        yAxis.setLabel("Controls"+(StringUtil.isBlank(this.controlFields)?"":"["+this.controlFields+"]"));
		        final ScatterChart<Number, Number>   chart =  new ScatterChart<>(xAxis,yAxis);
		        for(final XYChart.Series<Number,Number> series:partition.getSeries())
			        {
					chart.getData().add(series);
			        }
				String title="Case/Control";
				if(!args.isEmpty())
					{
					title= args.get(0);
					int slash=title.lastIndexOf("/");
					if(slash!=-1) title=title.substring(slash+1);
					if(title.endsWith(".vcf.gz")) title=title.substring(0, title.length()-7);
					if(title.endsWith(".vcf")) title=title.substring(0, title.length()-4);
					}
				if(userTitle!=null) title=userTitle;
				chart.setTitle(title);
				chart.setAnimated(false);
				//chart.setLegendSide(this.legendSide);
				
				
				final RExporter rExporter = new RExporter();
				final PrintWriter pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
				rExporter.exportToR(pw, chart);
				pw.flush();
				pw.close();

				
				}
			catch(final Exception err)
				{	
				LOG.error(err);
				return -1;
				}
			finally {
				CloserUtil.close(in);
				}
			
			
	        
	       
	        return 0;
	        }
		
	
	public static void main(final String[] args) {
		new CaseControlJfx().instanceMainWithExit(args);
	}
	
	
	
}
