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
package com.github.lindenb.jvarkit.tools.vcfgroupbypop;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.JVarkitVersion;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

@Program(name="vcfgroupbypop",
	description="create INFO data by population",
	creationDate="20190319",
	modificationDate="20230712",
	keywords= {"vcf","pedigree","population"}
	)
public class VcfGroupByPopulation extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfGroupByPopulation.class);
	@Parameter(names={"-p","--mapping","--sample2pop"},description="mapping file: each line is (SAMPLE)\\t(POP)\\n",required=true)
	private Path mappingFile=null;
	@Parameter(names={"-m","--min-fisher"},description="min inclusive value of fisher test. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double min_fisher = 0.0 ;
	@Parameter(names={"-M","--max-fisher"},description="max inclusive value of fisher test. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double max_fisher = 1.0 ;

	private final List<Pop> populations=new ArrayList<>();

	
	private static class Pop {
		final String name;
		final Set<String> samples = new HashSet<>();
		final VCFInfoHeaderLine ac;
		final VCFInfoHeaderLine an;
		final VCFInfoHeaderLine af;
		final VCFInfoHeaderLine dp;
		final VCFInfoHeaderLine n_missing;
		final VCFInfoHeaderLine n_HET;
		final VCFInfoHeaderLine n_HOM_REF;
		final VCFInfoHeaderLine n_HOM_VAR;
		
		Pop(final String name) {
			this.name=name;
			this.ac = new VCFInfoHeaderLine(name+"_"+VCFConstants.ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer,VCFConstants.ALLELE_COUNT_KEY+" for population "+name );
			this.an = new VCFInfoHeaderLine(name+"_"+VCFConstants.ALLELE_NUMBER_KEY, 1, VCFHeaderLineType.Integer,VCFConstants.ALLELE_NUMBER_KEY+" for population "+name );
			this.af = new VCFInfoHeaderLine(name+"_"+VCFConstants.ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer,VCFConstants.ALLELE_FREQUENCY_KEY+" for population "+name );
			this.dp = new VCFInfoHeaderLine(name+"_"+VCFConstants.DEPTH_KEY,1, VCFHeaderLineType.Integer,VCFConstants.DEPTH_KEY+" for population "+name );
			this.n_missing = new VCFInfoHeaderLine(name+"_MISSING",1, VCFHeaderLineType.Integer,"Number of missing genotypes for population"+name);
			this.n_HET = new VCFInfoHeaderLine(name+"_HET",1, VCFHeaderLineType.Integer,"Number of HET genotypes for population"+name);
			this.n_HOM_REF = new VCFInfoHeaderLine(name+"_HOM_REF",1, VCFHeaderLineType.Integer,"Number of HOM_REF genotypes for population"+name);
			this.n_HOM_VAR = new VCFInfoHeaderLine(name+"_HOM_VAR",1, VCFHeaderLineType.Integer,"Number of HOM_VAR genotypes for population"+name);
		}
		List<VCFInfoHeaderLine> getInfos() {
			return Arrays.asList(ac,an,af,dp,n_missing,n_HET,n_HOM_REF,n_HOM_VAR);
		}
		
		
		int[] fisher(final VariantContext ctx) {
			int R=0;
			int A=0;
			for(String sn:this.samples) {
				final Genotype g = ctx.getGenotype(sn);
				if(g.getAlleles().stream().allMatch(al->al.isReference() || al.isNoCall())){
					R++;
					}
				else
					{
					A++;
					}
				}
			return new int[] {R,A};
			}
		
		void annotate(final VariantContext ctx,VariantContextBuilder vcb) {
			final List<Allele> alts  = ctx.getAlternateAlleles();
			final List<Integer> ac_array = new ArrayList<>(alts.size());
			final List<Double> af_array = new ArrayList<>(alts.size());
			boolean dp_flag=false;
			int an=0;
			int dp = 0;
			int n_missing=0;
			int n_RR=0;
			int n_AR=0;
			int n_AA=0;
			for(String sn:this.samples) {
				final Genotype g = ctx.getGenotype(sn);
				switch(g.getType()) {
					case HOM_REF: n_RR++; break;
					case HET : n_AR++;break;
					case NO_CALL: n_missing++;break;
					case HOM_VAR: n_AA++;break;
					default: break;
					}
				for(Allele a: g.getAlleles()) {
					if(a.isNoCall()) continue;
					an++;
					}
				if(g.hasDP()) {
					dp_flag = true;
					dp += g.getDP();
					}
				}
			vcb.attribute(this.n_missing.getID(), n_missing);
			vcb.attribute(this.n_HOM_REF.getID(), n_RR);
			vcb.attribute(this.n_HET.getID(), n_AR);
			vcb.attribute(this.n_HOM_VAR.getID(), n_AA);

			vcb.attribute(this.an.getID(), an);
			if(dp_flag) vcb.attribute(this.dp.getID(), dp);
			
			for(int i=0;i< alts.size();i++)
				{
				int ac =0;
				for(String sn:this.samples) {
					final Genotype g = ctx.getGenotype(sn);
					for(Allele a: g.getAlleles()) {
						if(!a.equals(alts.get(i))) continue;
						ac++;
						}
					}
				ac_array.add(ac);
				if(an>0) af_array.add( ac/(double)an);
				}
			vcb.attribute(this.ac.getID(),ac_array);
			if(an>0) vcb.attribute(this.af.getID(),af_array);
			}
	}

	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator vcfIn, VariantContextWriter out) {
				try {
					final VCFHeader header= vcfIn.getHeader();
				
				final Set<String> samplesInVcf=new HashSet<>( header.getSampleNamesInOrder());
				for(Pop pop: this.populations) {
					pop.samples.retainAll(samplesInVcf);
					}
				this.populations.removeIf(P->P.samples.isEmpty());
				
				for(String sample: header.getSampleNamesInOrder())
					{
					if(this.populations.stream().noneMatch(P->P.samples.contains(sample)))
						{
						LOG.warn("Sample "+sample+" not affected to a population");
						}
					}
				
				final Set<VCFHeaderLine> metaData=new LinkedHashSet<>();
		
				for(Pop pop:this.populations) {
					metaData.addAll(pop.getInfos());
					}
				
				for(int i=0;i< populations.size();i++) {
						for(int j=i+1;j< populations.size();j++) {
						metaData.add( new VCFInfoHeaderLine(
								populations.get(i).name+"_"+populations.get(j).name+"_FISHER",1,VCFHeaderLineType.Float,"Fisher Test"
								));
						}
					}
				
				
				JVarkitVersion.getInstance().addMetaData(this, header);
				out.writeHeader(header);
						
				while(vcfIn.hasNext())
					{
					final VariantContext ctx= vcfIn.next();
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				
					for(Pop pop: this.populations) {
						pop.annotate(ctx, vcb);
						}
					
					
					for(int i=0;i +1 < populations.size();i++) {
						final Pop popi = populations.get(i);
						final int[] RAi =popi.fisher(ctx);
						
						for(int j=i+1;j< populations.size();j++) {
							final Pop popj = populations.get(j);
							final int[] RAj =popj.fisher(ctx);
							final FisherExactTest fisherTest = FisherExactTest.compute(
									RAi[1],RAi[0],
									RAj[1],RAj[0]
									);
							final double fisher = fisherTest.getAsDouble();
							if(fisher> this.max_fisher) continue;
							if(fisher< this.min_fisher) continue;
							vcb.attribute(popi.name+"_"+popj.name+"_FISHER",fisher);
							}
						}
		
					
					out.add(vcb.make());
					}
				return 0;
				}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected int beforeVcf()
		{
		try(BufferedReader br = IOUtils.openPathForBufferedReading(this.mappingFile)) {
			final Map<String,Pop> name2pop = new TreeMap<>();
			String line;
			while((line=br.readLine())!=null)
					{
					if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
					int space=line.indexOf('\t');
					if(space<=0) throw new IOException("tab missing in "+line);
					final String sample = line.substring(0,space);
					final String popName= line.substring(space+1).trim();
					if(StringUtils.isBlank(sample))  throw new IOException("empty sample in "+line);
					if(StringUtils.isBlank(popName))  throw new IOException("empty pop in "+line);
					Pop prevpo= name2pop.get(popName);
					if(prevpo==null) {
						prevpo=new Pop(popName);
						name2pop.put(popName,prevpo);
						}
					prevpo.samples.add(sample);
					}
			this.populations.addAll(name2pop.values());
			}
		catch(final IOException err ) {
			LOG.error(err);
			return -1;
			}
		
		return super.beforeVcf();
		}
	
	
	public static void main(final String[] args)
		{
		new VcfGroupByPopulation().instanceMainWithExit(args);
		}
}
