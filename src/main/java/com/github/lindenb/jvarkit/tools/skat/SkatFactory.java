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
package com.github.lindenb.jvarkit.tools.skat;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.tools.burden.MafCalculator;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.ProcessExecutor;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

@XmlRootElement(name="skat")
public class SkatFactory {
	private static final Logger LOG = Logger.build(SkatFactory.class).make();

	@Parameter(names={"--skat-adjusted"},description="SKAT adjusted")
	private boolean adjusted = false;
	@Parameter(names={"--skat-optimized"},description="SKAT optimized (SKATO)/ davies method.")
	private boolean optimal = false;
	@Parameter(names={"--skat-random-seed"},description="Rstats value for `set.seed`. -1 == use random")
	private int set_random_seed_value = -1;
	@Parameter(names={"--skat-accept-filtered"},description="accept variants FILTER-ed")
	private boolean acceptFILTERED = false;
	@Parameter(names={"--skat-num-retry"},description="compute n-times the p-value")
	private int n_retry = 1;

	private String RScript= "Rscript";
	
	public static interface SkatExecutor {
		
		/** a filter to be used in the VcfReader to discard variant as soon as possible */
		public Predicate<VariantContext> getUpstreamVariantFilter();
		
		public SkatResult execute(
			List<VariantContext> variants,
			Collection<Pedigree.Person> ped
			);
		};
	
		
		
		
	public static interface SkatResult extends Comparable<SkatResult>
		{
		public String getMessage();
		public double getPValue();
		public boolean isError();
		@Override
		public default int compareTo(final SkatResult o)
			{
			return Double.compare(this.getPValue(), o.getPValue());
			}
		}
	
	private static class ResultImpl implements SkatResult
		{
		final double pvalue;
		ResultImpl(double pvalue) { this.pvalue = pvalue;}
		@Override public String getMessage() { return "";}
		@Override public boolean isError() { return false;}
		@Override public double getPValue() { return this.pvalue;}
		@Override
		public String toString() {
			return String.valueOf(this.pvalue);
			}
		
		}
	
	private static class ResultError implements SkatResult
		{
		final String message;
		ResultError(final String message) { this.message = StringUtil.isBlank(message)?"ERROR":message;}
		@Override public String getMessage() { return this.message;}
		@Override public boolean isError() { return true;}
		@Override public double getPValue() { throw new IllegalStateException("cannot get p-value for an error : "+this.message);}
		@Override
		public String toString() { return "ERROR:"+this.message;}
		}

public void setAdjusted(boolean adjusted) {
	this.adjusted = adjusted;
	}
@XmlElement(name = "adjusted")
public boolean isAdjusted() {
	return adjusted;
}

public void setOptimal(boolean optimal) {
	this.optimal = optimal;
}
@XmlElement(name = "optimal")
public boolean isOptimal() {
	return optimal;
	}

public SkatExecutor build() {
	return new ExecutorImpl();
	}


private class ExecutorImpl implements SkatExecutor {

private final boolean	adjusted = SkatFactory.this.adjusted;
private final boolean	optimal = SkatFactory.this.optimal;
private final boolean acceptFILTERED = SkatFactory.this.acceptFILTERED;
private final int n_retry = SkatFactory.this.n_retry;
private final String RScript=  SkatFactory.this.RScript;
private final int set_random_seed_value = SkatFactory.this.set_random_seed_value;
private final File scriptFile;
private final File saveFile;

public ExecutorImpl() {
	try 
		{
		this.scriptFile = File.createTempFile("skat", ".R");
		//this.scriptFile.deleteOnExit();
		this.saveFile = File.createTempFile("skat", ".txt");
		//this.saveFile.deleteOnExit();		
		if(this.n_retry<1) throw new IllegalArgumentException("n_retry <1");
		}
	catch(final IOException err)
		{
		throw new RuntimeIOException(err);
		}
	}

@Override
public Predicate<VariantContext> getUpstreamVariantFilter() {
		return new Predicate<VariantContext>()
			{
			@Override
			public boolean test(final VariantContext V) {
				if(!acceptFILTERED && V.isFiltered()) return false;
				if(V.getNAlleles()!=2) return false;
				return true;
				}
			};
		}


private MafCalculator calculateMaf(final VariantContext ctx,final Collection<Pedigree.Person>  samples) {
	final Allele observed_alt = ctx.getAltAlleleWithHighestAlleleCount();
	final MafCalculator mafCalculator= new MafCalculator(observed_alt,ctx.getContig());
	for (final Pedigree.Person person : samples) {
		final Genotype genotype = ctx.getGenotype(person.getId());
		mafCalculator.add(genotype,  person.isMale());
		}
	return mafCalculator;
	}


private boolean isAdjusted() {
	return adjusted;
	}

private boolean isOptimal() {
	return optimal;
	}


private String getMethod() {
	return isOptimal()?"optimal":"davies";
	}

private String getKernel() {
	return "linear.weighted";
	}

	
@Override
public SkatFactory.SkatResult execute(
		List<VariantContext> variants,
		final Collection<Pedigree.Person> ped
		)
	{
	if(variants==null || variants.isEmpty()) return new ResultError("no variant");
	if(ped==null || ped.isEmpty()) return new ResultError("ped is empty");
	variants = variants.stream().
			filter(this.getUpstreamVariantFilter()).
			collect(Collectors.toList())
			;
	if(variants.isEmpty()) return new ResultError("no valid variants");
	final Set<String> samplesInVcf = variants.
			stream().
			flatMap(V->V.getGenotypes().stream()).
			map(G->G.getSampleName()).
			collect(Collectors.toSet()
			);
	
	final List<Pedigree.Person> samples = ped.
		stream().
		filter(P->(P.isAffected() || P.isUnaffected()) && samplesInVcf.contains(P.getId())).
		collect(Collectors.toList())
		;
	if(samples.isEmpty()) return new ResultError("no valid persons");
	
	variants =  variants.stream().
			filter(V->!calculateMaf(V,samples).isEmpty()).
			collect(Collectors.toList());
	if(variants.isEmpty()) return new ResultError("no variants with valid MAF");
	
	PrintWriter pw = null;
	try {
		pw = new PrintWriter(this.scriptFile);
		final String[] command = new String[]{this.RScript,this.scriptFile.getAbsolutePath()};
		
		
		pw.println("library(\"SKAT\")");
		
		if(this.set_random_seed_value!=-1)
			{
			pw.println("set.seed("+this.set_random_seed_value+")");
			}
		
		pw.print("phenotypes <- c(");
		pw.print(samples.stream().map(P->P.isUnaffected()?"0":"1").collect(Collectors.joining(",")));
		pw.println(")");
		
		pw.print("MAFs <- c(");
		pw.print(variants.stream().map(V->String.valueOf(calculateMaf(V, samples).getMaf())).collect(Collectors.joining(",")));
		pw.println(")");


		pw.print("genotypes <- c("); boolean first=true;
		for(final VariantContext ctx:variants)
			{
			for(final Pedigree.Person p: samples)
				{
				final Genotype genotype= ctx.getGenotype(p.getId());
				if(!first) pw.print(","); 
				first=false;
				if(genotype.isHomVar())
					{
					pw.print(2);
					}
				else if(genotype.isHet())
					{
					pw.print(1);
					}
				else
					{
					pw.print(0);
					}
				}
			}
		pw.println(")");
		pw.println("genot_matrix <- matrix(genotypes,nrow = "+samples.size()+", byrow = FALSE)");
		
		pw.println("obj2 <- SKAT_Null_Model(phenotypes~1, out_type=\"D\", Adjustment="+(isAdjusted()?"T":"F")+")");
		
		
		pw.println( "vector_values <- c()");
		pw.println( "n <- "+this.n_retry);
		pw.println( "while( n > 0 ) {");
		pw.println( "skat_value =SKAT(Z= genot_matrix ,weights=1/sqrt("+samples.size()+"*MAFs*(1-MAFs)),obj=obj2, kernel=\""+
				this.getKernel() + "\", method=\""+ this.getMethod()+"\")");
		pw.println( " vector_values <- c( vector_values ,c(skat_value$p.value )) ");
		pw.println( " n <- n - 1");
		pw.println( " }");
		

		pw.println("cat(median(vector_values),file=\""+ this.saveFile.getPath()+"\")");
		pw.println("quit()");
		pw.flush();
		pw.close();
		final int ret = ProcessExecutor.execute(command);
		if(ret!=0)
			{
			return new ResultError("RScript returned "+ret+" (!=0).");
			}
		try
			{
			return new ResultImpl(Double.parseDouble(IOUtil.slurp(this.saveFile)));
			}
		catch(final Throwable err2)
			{
			LOG.error(err2);
			return  new ResultError(err2.getMessage());
			}
		}
	catch(final Throwable err)
		{
		LOG.error(err);
		return  new ResultError(err.getMessage());
		}
	finally {
		CloserUtil.close(pw);
		this.scriptFile.delete();
		this.saveFile.delete();
		}
	}
}
}
