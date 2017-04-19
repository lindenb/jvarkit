/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

@Program(name="vcfrenamechr",description="Convert the names of the chromosomes in a VCF file")
public class ConvertVcfChromosomes extends com.github.lindenb.jvarkit.util.jcommander.Launcher {
	private static final Logger LOG = Logger.build(ConvertVcfChromosomes.class).make();
	
	
	@Parameter(names={"-c","-convert"},description="What should I do when  a converstion is not found")
	private ContigNameConverter.OnNotFound onNotFound=ContigNameConverter.OnNotFound.RAISE_EXCEPTION;
	@Parameter(names={"-f","--mapping","-m"},description="load a custom name mapping. Format (chrom-source\\tchrom-dest\\n)+",required=true)
	private File mappingFile=null;
	@Parameter(names={"-o","--out"},description="output vcf")
	private File output= null;

	
	
	public ConvertVcfChromosomes()
		{
		}
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator iterin, VariantContextWriter out) {
		final ContigNameConverter customMapping = ContigNameConverter.fromFile(mappingFile);
		customMapping.setOnNotFound(this.onNotFound);
		final Set<String> unseen=new HashSet<>();
		final VCFHeader header1=iterin.getHeader();
		final VCFHeader header2=new VCFHeader(
				header1.getMetaDataInInputOrder().stream().
					filter(L->!L.getKey().equals(VCFHeader.CONTIG_KEY)).collect(Collectors.toSet()),
					header1.getSampleNamesInOrder()
					);

		if(header1.getSequenceDictionary()!=null)
			{
			header2.setSequenceDictionary(customMapping.convertDictionary(header1.getSequenceDictionary()));
			}
		out.writeHeader(header2);
		while(iterin.hasNext()) {
			final VariantContext ctx=iterin.next();
			final String newName= customMapping.apply(ctx.getContig());
			if(newName==null)
				{
				if(unseen.size()<1000 && !unseen.contains(ctx.getContig()))
					{
					LOG.warn("Cannot find contig for "+ctx.getContig());
					unseen.add(ctx.getContig());
					}
				//skip unknown chromosomes
				continue;
				}
			final VariantContextBuilder vcb= new VariantContextBuilder(ctx);
			vcb.chr(newName);
			out.add(vcb.make());
			}
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		if( this.mappingFile==null)
			{
			throw new JvarkitException.CommandLineError("undefined mapping file");
			}
		
		return doVcfToVcf(args, this.output);
		}

	public static void main(String[] args)
		{
		new ConvertVcfChromosomes().instanceMainWithExit(args);
		}
	}
