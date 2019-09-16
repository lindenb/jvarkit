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
package com.github.lindenb.jvarkit.tools.structvar.manta;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.samtools.Decoy;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.sv.StructuralVariantComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;

/**
 BEGIN_DOC
 
 
 END_DOC

 */

@Program(name="mantafisher",
description="Fisher test from Manta Data",
keywords= {"sv","burden","manta"},
generate_doc=false,
creationDate="20190916",
modificationDate="20190916"
)
public class MantaFisher extends Launcher {
	private static final Logger LOG = Logger.build( MantaFisher.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private StructuralVariantComparator svComparator = new StructuralVariantComparator();
	@Parameter(names={"-c","--contig"},description="limit to this contig")
	private String limitContig = null;
	@Parameter(names={"--no-bnd"},description="discar BND")
	private boolean discard_bnd = false;

	
	private class SVKey {
		final VariantContext archetype;
		SVKey(final VariantContext archetype) {
			this.archetype = archetype;
			}
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof SVKey)) return false;
			final SVKey other = SVKey.class.cast(obj);
			return svComparator.test(this.archetype, other.archetype);
			}
		@Override
		public int hashCode() {
			int h=1;
			h= h*31 + archetype.getContig().hashCode();
			h= h*31 + archetype.getStart();
			if(svComparator.isTestingSvTypes()) {
				h= h*31 + this.archetype.getAttributeAsString(VCFConstants.SVTYPE,"").hashCode();
				}
			return h;
			}
		@Override
		public String toString() {
				return new SimpleInterval(this.archetype).toString();
			}
		}
	
	private static class SVValue {
		int count_affected_wild = 0;
		int count_affected_alt = 0;
		int count_unaffected_wild = 0;
		int count_unaffected_alt = 0;
		/** temporary flag, set if value was seen when scanning a VCF file */
		boolean visited_flag = false;
		}

	private static class VcfInput {
		Path vcfPath;
		boolean affected= false;
		int contigCount=0;
		}
	
@Override
public int doWork(final List<String> args) {
	PrintWriter out = null;
	try {
		final List<VcfInput> inputs = new ArrayList<>();
		SAMSequenceDictionary dict=null;
		final String input = oneFileOrNull(args);
		try(BufferedReader br=(input==null?
				IOUtils.openStreamForBufferedReader(stdin()):
				IOUtils.openURIForBufferedReading(oneFileOrNull(args))))
			{
			String line;
			while((line=br.readLine())!=null) {
				if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
				final String tokens[]= CharSplitter.TAB.split(line);
				if(tokens.length<2) throw new JvarkitException.TokenErrors(2, tokens);
				final VcfInput vcfInput = new VcfInput();
				vcfInput.vcfPath  = Paths.get(tokens[0]);
				IOUtil.assertFileIsReadable(vcfInput.vcfPath);
				final SAMSequenceDictionary dict1= SequenceDictionaryUtils.extractRequired(vcfInput.vcfPath);
				if(dict==null) {
					dict=dict1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict1))
					{
					throw new JvarkitException.DictionariesAreNotTheSame(dict1, dict);
					}
				if(tokens[1].equals("1")) {
					vcfInput.affected =true;
					}
				else if(tokens[1].equals("0")) {
					vcfInput.affected =false;
					}
				else
					{
					LOG.error("bad status in "+line);
					return -1;
					}
				inputs.add(vcfInput);
				}
			}
		if(inputs.isEmpty()) {
			LOG.error("no input found");
			return -1;
			}
		if(!StringUtils.isBlank(this.limitContig) && dict.getSequence(this.limitContig)==null) {
			LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(this.limitContig, dict));
			return -1;
			}
		if(inputs.stream().noneMatch(F->F.affected)) {
			LOG.error("No affected vcf");
			return -1;
			}
		if(inputs.stream().noneMatch(F->!F.affected)) {
			LOG.error("No non-affected vcf");
			return -1;
			}
		out = super.openPathOrStdoutAsPrintWriter(this.outputFile);
		out.print("#contig");
		out.print("\t");
		out.print("start-0");
		out.print("\t");
		out.print("end");
		out.print("\t");
		out.print("length");
		out.print("\t");
		out.print("svtype");
		out.print("\t");
		out.print("fisher");
		out.print("\t");
		out.print("case-ref");
		out.print("\t");
		out.print("case-alt");
		out.print("\t");
		out.print("control-ref");
		out.print("\t");
		out.print("control-alt");
		out.println();
		
		final Decoy decoy = Decoy.getDefaultInstance();
		for(final SAMSequenceRecord ssr: dict.getSequences()) {
			if(!StringUtils.isBlank(this.limitContig)) {
				if(!ssr.getSequenceName().equals(this.limitContig)) continue;
				}
			
			LOG.info("contig "+ssr.getSequenceName());
			if(decoy.isDecoy(ssr.getSequenceName())) continue;
			final Map<SVKey,SVValue> variants2count = new HashMap<>();
			for(final VcfInput vcfinput: inputs) {
				vcfinput.contigCount=0;//reset count for this contig
				
				try(VCFFileReader vcfFileReader = new VCFFileReader(vcfinput.vcfPath,true)) {
					vcfFileReader.query(ssr.getSequenceName(), 1, ssr.getSequenceLength()).
						stream().
						filter(V->discard_bnd==false || !V.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")).
						map(V->new VariantContextBuilder(V).
								unfiltered().
								noID().
								log10PError(VariantContext.NO_LOG10_PERROR).
								noGenotypes().
								rmAttribute("HOMSEQ").
								rmAttribute("SVINSSEQ").
								rmAttribute("LEFT_SVINSSEQ").
								rmAttribute("RIGHT_SVINSSEQ").
								make()).
						forEach(V->{
							variants2count.put(new SVKey(V), new SVValue());
							vcfinput.contigCount++;
						});
					
					LOG.info("contig "+ssr.getSequenceName()+" "+vcfinput.vcfPath+" N="+variants2count.size());
					}
				}
			if(variants2count.isEmpty()) continue;
			
			// build an interval tree for a faster access
			final IntervalTree<SVKey> intervalTree = new IntervalTree<>();
			for(final SVKey key: variants2count.keySet()) {
				final SimpleInterval r = new SimpleInterval(key.archetype).
						extend(this.svComparator.getBndDistance()+1);
				intervalTree.put(r.getStart(), r.getEnd(), key);
				}

			
			for(final VcfInput vcfinput: inputs) {
				LOG.info("now scanning "+ssr.getSequenceName()+" "+vcfinput.vcfPath+" containing "+ vcfinput.contigCount+" (total N="+variants2count.size()+")");
				try(VCFFileReader vcfFileReader = new VCFFileReader(vcfinput.vcfPath,true)) {
					//reset visited flag for this vcf
					variants2count.values().stream().forEach(V->V.visited_flag=false);
					
					final CloseableIterator<VariantContext> iter = vcfFileReader.query(ssr.getSequenceName(), 1, ssr.getSequenceLength());
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						
						if(this.discard_bnd && ctx.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")) continue;

						
						final Iterator<IntervalTree.Node<SVKey>> nodeIter = intervalTree.overlappers(ctx.getStart(),ctx.getEnd());
						while( nodeIter.hasNext()) {
							final SVKey key1 = nodeIter.next().getValue();
							if(!this.svComparator.test(key1.archetype, ctx)) continue;
  							final SVValue value = variants2count.get(key1);
							value.visited_flag = true;
							if(vcfinput.affected ) {
								value.count_affected_alt++;
								}
							else if(!vcfinput.affected) {
								value.count_unaffected_alt++;
								}
							}
						}
					iter.close();
					}
				// fill values that were not found here using 'visited_flag'
				for(final SVValue value : variants2count.values()) {
					if(value.visited_flag) continue;
					if(vcfinput.affected ) {
						value.count_affected_wild++;
						}
					else if(!vcfinput.affected) {
						value.count_unaffected_wild++;
						}
					}
				}
			
			for(final SVKey key : variants2count.keySet()) {
				final SVValue value = variants2count.get(key);
				out.print(key.archetype.getContig());
				out.print("\t");
				out.print(key.archetype.getStart()-1);
				out.print("\t");
				out.print(key.archetype.getEnd());
				out.print("\t");
				out.print(key.archetype.getLengthOnReference());
				out.print("\t");
				out.print(key.archetype.getAttributeAsString(VCFConstants.SVTYPE,"."));
				out.print("\t");
				out.print(value.count_affected_wild);
				out.print("\t");
				out.print(value.count_affected_alt);
				out.print("\t");
				out.print(value.count_unaffected_wild);
				out.print("\t");
				out.print(value.count_unaffected_alt);
				out.print("\t");
				out.print(FisherExactTest.compute(
						value.count_affected_wild, value.count_affected_alt,
						value.count_unaffected_wild, value.count_unaffected_alt).getAsDouble()
						);
				out.println();
				}
			}
		out.flush();
		out.close();
		out=null;
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(out);
		}
	}

public static void main(final String[] args) {
	new MantaFisher().instanceMainWithExit(args);
	}

}
