/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcf2r;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC



END_DOC
*/




@Program(name="vcf2r",
	description="Convert VCF to R so it can be used for burden testing",
	keywords={"vcf","burden","fisher"},
	modificationDate = "20240607",
	jvarkit_amalgamion = true
	)
public class VcfToRScript
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfToRScript.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;


	@Parameter(names={"-f","--function"},description="User defined R function to be called after each VCF")
	private String userDefinedFunName = "";

	
	@Parameter(names={"--info"},description="comma separated of extra INFO fields included in the variant dataframe")
	private String infoStr = "";


	@ParametersDelegate
	private CasesControls casesControls0 = new CasesControls();
	
	
	
	public VcfToRScript()
		{
		}
	
	
	
	private void process(final String filename,final VCFIterator in,PrintWriter pw,final Collection<String> infoFields0) {
		/* 
		 mail matile April 11: requires calling R for SKAT: needs
		 matrix of samples and genotypes */
		try {
			final VCFHeader header= in.getHeader();
			final CasesControls casesControls= this.casesControls0.clone();
			casesControls.retain(header);
			casesControls.checkHaveCasesControls();
			final List<VariantContext> variants = new ArrayList<>();
			final Set<String> all_samples = casesControls.getAll();
		

			
			boolean first=true;
			pw.println("# BEGIN  VCF ##########################################");
			pw.println("# Title ");
			
			pw.println("if(!exists(\"vcf.title\")) {\n    vcf.title <- \""+StringUtils.escapeC(filename)+"\"\n    }");
			
			
			first=true;
			pw.println("# samples ( 0: unaffected 1:affected)");
			
			pw.print("population <- data.frame(family=c(");	
			pw.print(all_samples.stream().map(S->StringUtils.doubleQuote(S)).collect(Collectors.joining(",")));
			pw.print("),name=c(");
			pw.print(all_samples.stream().map(S->StringUtils.doubleQuote(S)).collect(Collectors.joining(",")));
			pw.print("),status=c(");
			pw.print(all_samples.stream().map(S->casesControls.isControl(S)?"0":"1").collect(Collectors.joining(",")));
			pw.println("))");
						
			
			first=true;
			pw.println();
			pw.println("# genotypes as a list. Should be a multiple of length(samples).");
			pw.println("# 0 is homref (0/0), 1 is het (0/1), 2 is homvar (1/1)");
			pw.print("genotypes <- c(");
			while(in.hasNext()) {
				final VariantContext ctx =  in.next();
				if(!first) pw.print(",");
				pw.print(all_samples.stream().
					map(S->ctx.getGenotype(S)).
					map(genotype->{
						if (!genotype.hasAltAllele()) {
							return "0";
						} else if (genotype.hasRefAllele()) {
							return "1";
						} else {
							return "2";
						}}).
					collect(Collectors.joining(","))
					);
				first=false;
				variants.add(new VariantContextBuilder().noGenotypes().make());
				}// end reading vcf
			in.close();
			
			pw.println(")");
			first = true;
			
			pw.print("# variants.");
			pw.println();
			pw.print("variants <- data.frame(chrom=c(");
			pw.print(variants.stream().map(vc->StringUtils.doubleQuote(vc.getContig())).collect(Collectors.joining(",")));
			pw.print("),chromStart=c(");
			pw.print(variants.stream().map(vc->String.valueOf(vc.getStart())).collect(Collectors.joining(",")));
			pw.print("),chromEnd=c(");
			pw.print(variants.stream().map(vc->String.valueOf(vc.getEnd())).collect(Collectors.joining(",")));
			pw.print("),refAllele=c(");
			pw.print(variants.stream().map(vc->StringUtils.doubleQuote(vc.getReference().getDisplayString())).collect(Collectors.joining(",")));
			pw.print("),altAllele=c(");
			pw.print(variants.stream().map( vc->StringUtils.doubleQuote(vc.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(",")))).collect(Collectors.joining(",")));
			pw.print(")");
			
			for(String info: infoFields0) {
				final VCFInfoHeaderLine hdr = header.getInfoHeaderLine(info);
				if(hdr==null) continue;
				if(variants.stream().noneMatch(vc->vc.hasAttribute(info))) continue;
				pw.print(","+info+"=c(");
				if(hdr.getType()==VCFHeaderLineType.String) {
					pw.print(variants.stream().map( vc->StringUtils.doubleQuote(vc.getAttributeAsStringList(info,"").stream().collect(Collectors.joining(",")))).collect(Collectors.joining(",")));
					}
				else if(hdr.getType()==VCFHeaderLineType.Flag) {
					pw.print(variants.stream().map( vc->vc.hasAttribute(info)?"1":"0").collect(Collectors.joining(",")));
					}
				else if(hdr.getType()==VCFHeaderLineType.Float) {
					pw.print(variants.stream().map( vc->String.valueOf(vc.getAttributeAsDouble(info, 0.0))).collect(Collectors.joining(",")));
					}
				else if(hdr.getType()==VCFHeaderLineType.Integer) {
					pw.print(variants.stream().map( vc->String.valueOf(vc.getAttributeAsInt(info, 1))).collect(Collectors.joining(",")));
					}
				else
					{
					throw new IllegalStateException(info);
					}
				pw.print(")");
				}
			
			
			pw.println(")");
			
			if(!variants.isEmpty())
				{
				pw.println("# assert sizes");
				pw.println("stopifnot( length(genotypes) %% NROW(population) == 0 )");
				pw.println("stopifnot(NROW(variants) * NROW(population) == length(genotypes) )");
				
				
				if(this.userDefinedFunName==null || this.userDefinedFunName.trim().isEmpty()) {
					pw.println("## WARNING not user-defined R function was defined");
					}
				else
					{
					pw.println("# consumme data with user-defined R function ");
					pw.println(this.userDefinedFunName+"()");
					}
				}
			else
				{
				LOG.warn("No Variant found");
				}
			
			
			pw.println("# END VCF ##########################################");
			}
		catch(Throwable err) {
			throw new RuntimeException(err);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			
			this.casesControls0.load();
			this.casesControls0.checkHaveCasesControls();
			
			final List<Path> paths = IOUtils.unrollPaths(args);

			final List<String> infoFiels = Arrays.stream(this.infoStr.split("[ ,;\t]")).
					filter(it->!StringUtils.isBlank(it)).
					collect(Collectors.toSet()).
					stream().
					collect(Collectors.toList());
			
			
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				pw.println("# This is the header of the generated R script. ");
				pw.println("# generated by "+ getProgramName());
				pw.println("# version "+getGitHash());
				
				if(!StringUtils.isBlank(this.userDefinedFunName)) {
					pw.println("# user defined function "+userDefinedFunName+" should have been inserted BEFORE this header. ");
					pw.println("#  something like 'cat user_function.R this_file.R |  R --no-save > result.txt'");
					}
				if(paths.isEmpty()) {
						try(VCFIterator r=VCFUtils.createVCFIteratorStdin()) {
							process("<STDIN>",r,pw,infoFiels);
							}
						}
					else
						{
						for(Path path:paths) {
							try(VCFIterator r=VCFUtils.createVCFIteratorFromPath(path)) {
								process(path.toString(),r,pw,infoFiels);
								}
							}
						}
				pw.flush();
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	public static void main(String[] args)
		{
		new VcfToRScript().instanceMainWithExit(args);
		}
	}
