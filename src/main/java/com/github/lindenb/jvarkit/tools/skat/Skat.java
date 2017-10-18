package com.github.lindenb.jvarkit.tools.skat;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.github.lindenb.jvarkit.tools.burden.MafCalculator;
import com.github.lindenb.jvarkit.util.Pedigree;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.ProcessExecutor;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

@XmlRootElement(name="skat")
public class Skat {
	private boolean adjusted = false;
	private boolean optimal = false;
	private boolean acceptFILTERED = false;
	private String RScript= "Rscript";
	
	public static interface SkatResult
		{
		public String getMessage();
		public double getPValue();
		public boolean isError();
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
		@Override public double getPValue() { throw new IllegalStateException("cannot get p-value for an error");}
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

private String getMethod() {
	return isOptimal()?"optimal":"davies";
	}

private String getKernel() {
	return "linear.weighted";
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

public Skat.SkatResult execute(
		List<VariantContext> variants,
		final Pedigree ped
		)
	{
	if(variants==null || variants.isEmpty()) return new ResultError("no variant");
	if(ped.isEmpty()) return new ResultError("ped is empty");
	variants = variants.stream().
			filter(V->this.acceptFILTERED || !V.isFiltered()).
			filter(V->V.getNAlleles()==2).
			collect(Collectors.toList());
	if(variants.isEmpty()) return new ResultError("no valid variants");
	final Set<String> samplesInVcf = variants.
			stream().
			flatMap(V->V.getGenotypes().stream()).
			map(G->G.getSampleName()).
			collect(Collectors.toSet()
			);
	
	final List<Pedigree.Person> samples = ped.getPersons().
		stream().
		filter(P->(P.isAffected() || P.isUnaffected()) && samplesInVcf.contains(P.getId())).
		collect(Collectors.toList())
		;
	if(samples.isEmpty()) return new ResultError("no valid persons");
	
	variants =  variants.stream().
			filter(V->!calculateMaf(V,samples).isEmpty()).
			collect(Collectors.toList());
	if(variants.isEmpty()) return new ResultError("no variants with valid MAF");
	
	File scriptFile = null;
	File saveFile = null;
	PrintWriter pw = null;
	try {
		
		scriptFile = File.createTempFile("skat", ".R");
		scriptFile.deleteOnExit();
		saveFile = File.createTempFile("skat", ".txt");
		saveFile.deleteOnExit();
		
		pw = new PrintWriter(scriptFile);
		final String[] command = new String[]{this.RScript,scriptFile.getAbsolutePath()};
		
		
		pw.println("library(\"SKAT\")");
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
		
		
		pw.println( "skat_value =SKAT(Z= genot_matrix ,weights=1/sqrt("+samples.size()+"*MAFs*(1-MAFs)),obj=obj2, kernel=\""+
				this.getKernel() + "\", method=\""+ this.getMethod()+"\")");
		
		pw.println("cat(skat_value$p.value,file=\""+ saveFile.getPath()+"\")");
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
			return new ResultImpl(Double.parseDouble(IOUtil.slurp(saveFile)));
			}
		catch(final Throwable err2)
			{
			err2.printStackTrace();
			return  new ResultError(err2.getMessage());
			}
		}
	catch(final Throwable err)
		{
		err.printStackTrace();
		return  new ResultError(err.getMessage());
		}
	finally {
		CloserUtil.close(pw);
		if(scriptFile!=null) scriptFile.delete();
		if(saveFile!=null) saveFile.delete();
	
		}
	}

public static void main(String[] args) {
	try {
		VCFFileReader r= new VCFFileReader(new File(args[0]),false);
		Pedigree ped = new Pedigree.Parser().parse(r.getFileHeader());
		List<VariantContext> variants= new ArrayList<>();
		CloseableIterator<VariantContext> iter=r.iterator();
		while(iter.hasNext()) variants.add(iter.next());
		iter.close();
		r.close();
		Skat skat=new Skat();
		SkatResult rez=skat.execute(variants, ped);
		System.err.println(rez);
		System.exit(0);
	} catch(Exception err) {
		err.printStackTrace();
	}
}
}
