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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.sv.StructuralVariantComparator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC
 
# Motivation
 
Find SV present in cases and controls.


 
# Input
 
 input is a tab delimited file with a header:
 
  - vcf: the path to the vcf
  - sample : the sample name associated to the vcf. If empty the name will be searched in the CHROM header
  - status : case OR control
 
 
# Example
 
 ```

 $ java -jar jvarkit.jar svcasescontrols input.tsv > output.vcf
 ```

END_DOC

 */

@Program(name="svcasescontrols",
description="Find SV present in cases but not in controls.",
keywords= {"sv","manta","vcf"},
creationDate="20240513",
modificationDate="20250901",
jvarkit_amalgamion = true,
menu="VCF Manipulation"
)
public class SVCasesControls extends Launcher {
	private static final Logger LOG = Logger.of( SVCasesControls.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private StructuralVariantComparator svComparator = new StructuralVariantComparator();
	@Parameter(names={"-c","--contig","--chrom"},description="limit to this contig")
	private String limitContig = null;
	@Parameter(names={"--no-bnd"},description="discar BND")
	private boolean discard_bnd = false;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();
	@Parameter(names={"--debug"},hidden=true)
	boolean debug_flag =false;
	@Parameter(names={"--complex"},description="By default this tool select the SV that are found in the cases but not in the controls. When using this flag, all variants for cases are extracted and a count of CASE/CONTROL having the SV is added in the INFO column.")
	private boolean complex_analysis=false;

	
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

	private String toString(final VariantContext VC) {
		return VC.getContig()+":"+VC.getStart()+"-"+VC.getEnd()+":"+VC.getLengthOnReference()+":"+VC.getAttributeAsString(VCFConstants.SVTYPE, ".");
		}

	private VCFReader open(final Path path)  {
		return VCFReaderFactory.makeDefault().open(path, true);
		}
	private Locatable extend(final VariantContext ctx) {
		if("BND".equals(ctx.getAttributeAsString(VCFConstants.SVTYPE,""))) {
			int x=1+svComparator.getBndDistance();
			return new SimpleInterval(ctx.getContig(),Math.max(1, ctx.getStart()+x),ctx.getEnd()+x);
			}
		return ctx;
		}

	@Override
	public int doWork(final List<String> args) {
		try {
			SAMSequenceDictionary dict=null;
			final List<VcfInput> vcfInputs = new ArrayList<>();

			final String input =  oneFileOrNull(args);
			try(BufferedReader br = super.openBufferedReader(input)) {
				String line = br.readLine();
				if(line==null) throw new IOException("cannot read first line of input "+(input==null?"stdin":input) );
				final FileHeader fileHeader = new FileHeader(line,S->Arrays.asList(CharSplitter.TAB.split(S)));
				fileHeader.assertColumnExists("vcf").assertColumnExists("status");
				while((line=br.readLine())!=null) {
					final FileHeader.RowMap row= fileHeader.toMap(line);
					final Path path = Paths.get(row.get("vcf"));
					String sample = row.getOrDefault("sample", "");
					String status = row.get("status");
					final boolean is_case;
					if(status.equals("case")) {
						is_case=true;
						}
					else if(status.equals("control") || status.equals(".") ||StringUtils.isBlank(status) ) {
						is_case=false;
						}
					else
						{
						LOG.error("bad status in "+row);
						return -1;
						}
					
					try(VCFReader r = open(path)) {
						final VCFHeader header = r.getHeader();
						if(StringUtils.isBlank(sample)) {
							final List<String> samples = header.getGenotypeSamples();
							if(samples.isEmpty()) {
								LOG.error("no genotype in "+path+" using path as sample name");
								sample = path.toString();
								}
							else if(samples.size()!=1) {
								LOG.error("expected one and only one sample in "+path+" but got "+samples.size());
								return -1;
								}
							else
								{
								sample = samples.get(0);
								}
							}
						final SAMSequenceDictionary dict1= SequenceDictionaryUtils.extractRequired(header);
						if(dict==null) {
							dict=dict1;
							}
						else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict1))
							{
							throw new JvarkitException.DictionariesAreNotTheSame(dict1, dict);
							}
						}
					final String final_sn = sample;
					if(vcfInputs.stream().anyMatch(V->V.sample.equals(final_sn))) {
						LOG.warning("DUPLICATE SAMPLE for sample "+final_sn);
						}
					vcfInputs.add(new VcfInput(path,sample,is_case));
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
			if(complex_analysis) {
				return runComplex(dict,vcfInputs);
				}
			else
				{
				return runSimple(dict,vcfInputs);
				}
			}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	}
		
		
	private int runSimple(final SAMSequenceDictionary dict,final List<VcfInput> vcfInputs) throws IOException {
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
				try(VCFReader r = open(input.vcfPath)) {
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
							try(CloseableIterator<VariantContext> iter = r.query(extend(ctx))) {
								while(iter.hasNext()) {
									final VariantContext ctx2  = iter.next();
									if(discard_bnd && "BND".equals(ctx2.getAttributeAsString(VCFConstants.SVTYPE,""))) {
										continue;
										}
									if(this.svComparator.test(ctx,ctx2)) {
										found_same=true;
										if(debug_flag) {
											LOG.debug("SAME: ");
											LOG.debug(toString(ctx));
											LOG.debug(toString(ctx2));
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
	
private static class VariantCount {
	VariantContext vc;
	final List<String> cases_alt=new ArrayList<>();
	int case_ref=0;
	final List<String> ctrls_alt=new ArrayList<>();
	int ctrl_ref=0;
	VariantCount(final VariantContext vc) {
		this.vc=vc;
		}
	}
	
private int runComplex(final SAMSequenceDictionary dict,final List<VcfInput> vcfInputs) throws IOException {
	VariantContextWriter out=null;
	final VCFInfoHeaderLine infoNCaseAlt= new VCFInfoHeaderLine("N_CASE_ALT", 1, VCFHeaderLineType.Integer, "count case Samples with ALT ");
	final VCFInfoHeaderLine infoNCaseRef= new VCFInfoHeaderLine("N_CASE_REF", 1, VCFHeaderLineType.Integer, "count case Samples without ALT");
	final VCFInfoHeaderLine infoNCtrlAlt= new VCFInfoHeaderLine("N_CTRL_ALT", 1, VCFHeaderLineType.Integer, "count control Samples with ALT ");
	final VCFInfoHeaderLine infoNCtrlRef= new VCFInfoHeaderLine("N_CTRL_REF", 1, VCFHeaderLineType.Integer, "count control Samples without ALT");
	final VCFInfoHeaderLine infoAFCtrls = new VCFInfoHeaderLine("AF_CTRL", 1, VCFHeaderLineType.Float, "Frequency in controls");
	final VCFInfoHeaderLine infoAFCases = new VCFInfoHeaderLine("AF_CASE", 1, VCFHeaderLineType.Float, "Frequency in cases");
	final VCFInfoHeaderLine infoFisher = new VCFInfoHeaderLine("FISHER", 1, VCFHeaderLineType.Float, "Fisher ExactTest");
	final VCFInfoHeaderLine infoCasesSamples = new VCFInfoHeaderLine("CASES_SAMPLES", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "cases samples with ALT");
	final VCFInfoHeaderLine infoCtrlsSamples = new VCFInfoHeaderLine("CTRLS_SAMPLES", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "ctrl samples with ALT");
	final int count_cases = (int)vcfInputs.stream().filter(V->V.is_case).count();
	final int count_ctrls = (int)vcfInputs.stream().filter(V->!V.is_case).count();
	for(final SAMSequenceRecord ssr: dict.getSequences()) {
		if(!StringUtils.isBlank(this.limitContig)) {
			if(!ssr.getSequenceName().equals(this.limitContig)) continue;
			}
		final IntervalTreeMap<List<VariantContext>> treeMap=new IntervalTreeMap<>();
		for(Path path: vcfInputs.stream().filter(V->V.is_case).map(V->V.vcfPath).collect(Collectors.toList())) {
			try(VCFReader r = open(path)) {
				try(CloseableIterator<VariantContext> iter = r.query(ssr)) {
					while(iter.hasNext()) {
						final VariantContext ctx  = new VariantContextBuilder(iter.next()).noGenotypes().make();
						if(discard_bnd && "BND".equals(ctx.getAttributeAsString(VCFConstants.SVTYPE,""))) {
							continue;
							}

						final Interval rgn = new Interval(ctx);
						List<VariantContext> L= treeMap.get(rgn);
						if(L==null) {
							L=new ArrayList<>();
							treeMap.put(rgn,L);
							}
						if(L.stream().noneMatch(V->this.svComparator.test(V,ctx))) {
							L.add(ctx);
							}
						}
					}
				}
			}
		
		final List<VariantCount> all_intervals =treeMap.values().stream().flatMap(L->L.stream()).sorted(new ContigDictComparator(dict).createLocatableComparator()).map(V->new VariantCount(V)).collect(Collectors.toList());
		for(int i=0;i< vcfInputs.size() /*&&  !all_intervals.isEmpty() no, must create VCF output */;i++) {
			final VcfInput input = vcfInputs.get(i);
			LOG.info(""+(i+1)+"/"+vcfInputs.size()+" contig: "+ssr.getSequenceName()+" "+input.vcfPath+" "+ input.sample+" "+(input.is_case?"CASE":"CONTROL")+" ");
			try(VCFReader r = open(input.vcfPath)) {
				final VCFHeader header = r.getHeader();
				for(VariantCount ctx:all_intervals) {
					boolean found_same=false;
					try(CloseableIterator<VariantContext> iter = r.query(extend(ctx.vc))) {
						while(iter.hasNext()) {
							final VariantContext ctx2  = iter.next();
							if(discard_bnd && "BND".equals(ctx2.getAttributeAsString(VCFConstants.SVTYPE,""))) {
								continue;
								}
							if(this.svComparator.test(ctx.vc,ctx2)) {
								found_same=true;
								break;
								}
							}
						}
					if(found_same==true) {
						if(input.is_case) {
							ctx.cases_alt.add(input.sample);
							}
						else
							{
							ctx.ctrls_alt.add(input.sample);
							}
						}
					else
						{
						if(input.is_case) {
							ctx.case_ref++;
							}
						else
							{
							ctx.ctrl_ref++;
							}
						}
					}
				
				if(out==null) {
					final Set<VCFHeaderLine> meta = new LinkedHashSet<>(header.getMetaDataInInputOrder());
					meta.add(infoNCaseAlt);
					meta.add(infoNCaseRef);
					meta.add(infoNCtrlAlt);
					meta.add(infoNCtrlRef);
					meta.add(infoAFCases);
					meta.add(infoAFCtrls);
					meta.add(infoFisher);
					meta.add(infoCasesSamples);
					meta.add(infoCtrlsSamples);
					final VCFHeader header2=new VCFHeader(meta,Collections.emptyList());
					JVarkitVersion.getInstance().addMetaData(this, header);
					
					out = this.writingVariants.dictionary(header).open(this.outputFile);
					out.writeHeader(header2);
					}
				if(all_intervals.isEmpty()) break;
				} // end open VCF
			} // end all inputs
		
		int x=0;
		while(x < all_intervals.size()) {
			final VariantCount vc1 = all_intervals.get(x);
			if(x+1 < all_intervals.size() &&
					this.svComparator.test(vc1.vc,all_intervals.get(x+1).vc) &&
					vc1.cases_alt.equals(all_intervals.get(x+1).cases_alt) &&
					vc1.case_ref == all_intervals.get(x+1).case_ref &&
					vc1.ctrls_alt.equals(all_intervals.get(x+1).ctrls_alt) &&
					vc1.ctrl_ref == all_intervals.get(x+1).ctrl_ref
					) {
				all_intervals.remove(x+1);
				}
			else
				{
				++x;
				}
			}
		
		for(VariantCount ctx:all_intervals) {
			final VariantContextBuilder vcb  = new VariantContextBuilder(ctx.vc).
				attribute(infoNCaseAlt.getID(),ctx.cases_alt.size()).
				attribute(infoNCaseRef.getID(),ctx.case_ref).
				attribute(infoNCtrlAlt.getID(),ctx.ctrls_alt.size()).
				attribute(infoNCtrlRef.getID(),ctx.ctrl_ref).
				noGenotypes()
				;
			vcb.attribute(infoFisher.getID(), FisherExactTest.compute(
					ctx.cases_alt.size(),ctx.case_ref,
					ctx.ctrls_alt.size(),ctx.ctrl_ref
					).getAsDouble());

			
			if(count_cases>0) {
				vcb.attribute(infoAFCases.getID(), ctx.cases_alt.size()/(double)count_cases);
				vcb.attribute(infoCasesSamples.getID(), new ArrayList<>(new TreeSet<>( ctx.cases_alt)));
				}
			if(count_ctrls>0) {
				vcb.attribute(infoAFCtrls.getID(), ctx.ctrls_alt.size()/(double)count_ctrls);
				vcb.attribute(infoCtrlsSamples.getID(), new ArrayList<>(new TreeSet<>( ctx.ctrls_alt)));
				}
			out.add(vcb.make());
			}
		}// end contigs
	out.close();
	return 0;
	}

public static void main(final String[] args) {
	new SVCasesControls().instanceMainWithExit(args);
	}

}
