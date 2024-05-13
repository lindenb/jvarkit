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
package com.github.lindenb.jvarkit.tools.structvar;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.sv.StructuralVariantComparator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 BEGIN_DOC
 
 # Input
 
 input is a list of indexed vcf files or one file with the '.list' suffix containing the path to the vcfs
 
 
 # Example
 
 ```
 $ find src -name "manta*z" > jeter.list
 $ java -jar jvarkit.jar svcasescontrols jeter.list 2> /dev/null
 
 (...)
 
 
 
 ```
 END_DOC

 */

@Program(name="svcasescontrols",
description="Find SV present in cases but not in controls.",
keywords= {"sv","manta","vcf"},
creationDate="20240513",
modificationDate="20240513",
jvarkit_amalgamion = true,
menu="VCF Manipulation"
)
public class SVCasesControls extends Launcher {
	private static final Logger LOG = Logger.build( SVCasesControls.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private StructuralVariantComparator svComparator = new StructuralVariantComparator();
	@Parameter(names={"-c","--contig"},description="limit to this contig")
	private String limitContig = null;
	@Parameter(names={"--no-bnd"},description="discar BND")
	private boolean discard_bnd = false;
	@Parameter(names={"-cases","--cases"},description="samples's name for cases")
	private List<String> casesStr=new ArrayList<>();
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();
	@Parameter(names={"--debug"},hidden=true)
	boolean debug_flag =false;
	
	private static class VcfInput {
		final Path vcfPath;
		final String sample;
		final boolean is_case;
		VcfInput(final Path vcfPath,final String sample,boolean is_case) {
			this.vcfPath=vcfPath;
			this.sample=sample;
			this.is_case=is_case;
			}
		}
	
@Override
public int doWork(final List<String> args) {
	try {
		SAMSequenceDictionary dict=null;
		final Set<String> cases = this.casesStr.stream().flatMap(S->Arrays.stream(S.split("[ ,\t;]"))).filter(S->!StringUtils.isBlank(S)).collect(Collectors.toSet());
		final List<VcfInput> vcfInputs = new ArrayList<>();
		
		
		final Function<VariantContext,String> vc2str = VC->VC.getContig()+":"+VC.getStart()+"-"+VC.getEnd()+":"+VC.getLengthOnReference()+":"+VC.getAttributeAsString(VCFConstants.SVTYPE, ".");
		
		for(Path path: IOUtils.unrollPaths(args)) {
			try(VCFFileReader r = new VCFFileReader(path, true)) {
				final VCFHeader header = r.getHeader();
				final List<String> samples = header.getGenotypeSamples();
				if(samples.size()!=1) {
					LOG.error("expected one and only one sample in "+path+" but got "+samples.size());
					return -1;
					}
				final String sn = samples.get(0);
				final SAMSequenceDictionary dict1= SequenceDictionaryUtils.extractRequired(header);
				if(dict==null) {
					dict=dict1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict1))
					{
					throw new JvarkitException.DictionariesAreNotTheSame(dict1, dict);
					}
				if(vcfInputs.stream().anyMatch(V->V.sample.equals(sn))) {
					LOG.warning("DUPLICATE SAMPLE for sample "+sn);
					}
				vcfInputs.add(new VcfInput(path, sn,cases.isEmpty() || cases.contains(sn)));
				}
			}
		if(vcfInputs.isEmpty()) {
			LOG.error("no input");
			return -1;
			}
		if(vcfInputs.stream().noneMatch(V->V.is_case)) {
			LOG.error("no input for CASE");
			return -1;
			}
		
		
		
		Collections.sort(vcfInputs,(A,B)->Integer.compare(A.is_case?0:1, B.is_case?0:1));
		
		
		if(!StringUtils.isBlank(this.limitContig) && dict.getSequence(this.limitContig)==null) {
			LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(this.limitContig, dict));
			return -1;
			}
		
		VariantContextWriter out=null;
		for(final SAMSequenceRecord ssr: dict.getSequences()) {
			final List<VariantContext> remains = new ArrayList<>();
			if(!StringUtils.isBlank(this.limitContig)) {
				if(!ssr.getSequenceName().equals(this.limitContig)) continue;
				}
			for(int i=0;i< vcfInputs.size();i++) {
				final VcfInput input = vcfInputs.get(i);
				LOG.info(""+(i+1)+"/"+vcfInputs.size()+" contig: "+ssr.getSequenceName()+" "+input.vcfPath+" "+ input.sample+" "+(input.is_case?"CASE":"CONTROL")+" "+(i>0?"remains: "+remains.size():"."));
				if(i==0 && !input.is_case) throw new IllegalStateException();
				if(i>0 && remains.isEmpty()) break;
				try(VCFFileReader r=new VCFFileReader(input.vcfPath,true)) {
					final VCFHeader header = r.getHeader();
					if(i==0) {
						 try(CloseableIterator<VariantContext> iter = r.query(ssr)) {
							while(iter.hasNext()) {
								final VariantContext ctx  = iter.next();
								if(discard_bnd && "BND".equals(ctx.getAttributeAsString(VCFConstants.SVTYPE,""))) {
									continue;
									}
								remains.add(ctx);
								}
						 	}
						 if(out==null) {
							final VCFHeader header2=new VCFHeader(header.getMetaDataInInputOrder(),Collections.emptyList());
							JVarkitVersion.getInstance().addMetaData(this, header);
							
							out = this.writingVariants.dictionary(header).open(this.outputFile);
							out.writeHeader(header2);
							}
						}
					else
						{
						int x=0;
						while(x< remains.size()) {
							final VariantContext ctx  = remains.get(x);
							boolean found_same=false;
							try(CloseableIterator<VariantContext> iter = r.query(ctx)) {
								while(iter.hasNext()) {
									final VariantContext ctx2  = iter.next();
									if(discard_bnd && "BND".equals(ctx2.getAttributeAsString(VCFConstants.SVTYPE,""))) {
										continue;
										}
									if(this.svComparator.test(ctx,ctx2)) {
										found_same=true;
										if(debug_flag) {
											LOG.debug("SAME: ");
											LOG.debug(vc2str.apply(ctx));
											LOG.debug(vc2str.apply(ctx2));
											}
										break;
										}
									}
							 	}
							boolean remove_variant=(input.is_case && !found_same) || (!input.is_case && found_same);
							if(remove_variant) {
								remains.remove(x);
								}
							else
								{
								x++;
								}
							}
						}
					}
				}
			for(VariantContext ctx:remains) {
				out.add(new VariantContextBuilder(ctx).noGenotypes().make());
				}
			}
		out.close();
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	}

public static void main(final String[] args) {
	new SVCasesControls().instanceMainWithExit(args);
	}

}
