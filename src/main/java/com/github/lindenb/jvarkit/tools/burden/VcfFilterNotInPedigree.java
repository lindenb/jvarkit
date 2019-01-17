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


History:

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

## Example



END_DOC
 */
@Program(name="vcffilternotinpedigree",
	description="Adds a FILTER 'NotInPedigree' if the only not(homref) genotypes are not in a pedigree",
	keywords={"burden","vcf","pedigree"}
	)
public class VcfFilterNotInPedigree
	extends Launcher
	{

	private static final Logger LOG = Logger.build(VcfFilterNotInPedigree.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	
	@XmlType(name="vcffilternotinpedigree")
	@XmlRootElement(name="vcffilternotinpedigree")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{
			@XmlElement(name="filter")
			@Parameter(names={"-f","--filter"},description="FILTER name. Will be set for variant where the only genotypes non-homref are NOT in the pedigree")
			private String filterName = "NoGenotypeInPedigree";
		
			@XmlElement(name="discard")
			@Parameter(names={"-r","--remove"},description="remove the variant instead of setting the FILTER (hard filtering)")
			private boolean dicardVariant = false;
		
			@XmlElement(name="singleton")
			@Parameter(names={"-s","--singleton"},description="Variant is flagged/FILTERed as SingletonAlt if the ALT if found in less or equal times 'singleton-times' in the genotypes. -1:ignore")
			private int singleton = 1 ;
			
			@XmlElement(name="singleton-filter")
			@Parameter(names={"-sf","--sfilter"},description="FILTER name for option --singleton")
			private String singletonfilterName = "SingletonAlt";
			
			@XmlElement(name="pedigree")
			@Parameter(names={"-p","--pedigree"},description="[20180406]" + Pedigree.OPT_DESCRIPTION+" Default is to try to read the pedigree in the VCF header")
			private File pedigreeFile = null;
			@XmlElement(name="pedigree")
			@Parameter(names={"-gtf","--ignore-filtered-gt"},description="[20180406] Do not consider a *genotype* if it is FILTERED.")
			private boolean ignoreFilteredGT = false;

			
			private class CtxWriter extends DelegateVariantContextWriter
				{
				private VCFFilterHeaderLine singletonFilter = null;
				private VCFFilterHeaderLine filter = null;
				private Set<Pedigree.Person> individuals = null;
				private final int IGNORE_SINGLETON=-1;
				private final Predicate<Genotype> acceptFilteredGenotype = G -> 
					G!=null && (CtxWriterFactory.this.ignoreFilteredGT == false ||  !G.isFiltered())
					;
				
				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				@Override
				public void writeHeader(final VCFHeader header) {
					
					final Pedigree pedigree;
							
					if(CtxWriterFactory.this.pedigreeFile==null) {		
						pedigree = Pedigree.newParser().parse(header);
						}
					else
						{
						try {
							pedigree = Pedigree.newParser().parse(CtxWriterFactory.this.pedigreeFile);
						} catch (final IOException err) {
							throw new RuntimeIOException("Cannot read pedigree in file: "+CtxWriterFactory.this.pedigreeFile,err);
							}
						}
					if(pedigree.isEmpty()) {
						throw new JvarkitException.PedigreeError("No pedigree found in header. use VcfInjectPedigree to add it");
						}
					if(!pedigree.verifyPersonsHaveUniqueNames()) {
						throw new JvarkitException.PedigreeError("I can't use this pedigree in VCF because two samples have the same ID");
					}

					final Set<String> samplesNames= new HashSet<>(header.getSampleNamesInOrder());
					this.individuals = new HashSet<>(pedigree.getPersons());
					final Iterator<Pedigree.Person> iter= this.individuals.iterator();
					while(iter.hasNext())
					{
						final Pedigree.Person person = iter.next();
						if(!(samplesNames.contains(person.getId()))) {
							LOG.warn("Ignoring "+person+" because not in VCF header.");
							iter.remove();
						}
					}
					
					this.filter = new VCFFilterHeaderLine(
							CtxWriterFactory.this.filterName,
							"Will be set for variant where the only genotypes non-homref are NOT in the pedigree "
							);
					this.singletonFilter = new VCFFilterHeaderLine(
							CtxWriterFactory.this.singletonfilterName,
							"The ALT allele is found in less or equals than "+CtxWriterFactory.this.singleton+" individuals in the cases/controls"
							);
					
					final VCFHeader h2= new VCFHeader(header);
					h2.addMetaDataLine(filter);
					if(CtxWriterFactory.this.singleton!=IGNORE_SINGLETON) {
						h2.addMetaDataLine(this.singletonFilter);
						}					
					super.writeHeader(h2);
					}
				
				@Override
				public void add(final VariantContext ctx) {
					final boolean in_pedigree= this.individuals.stream().
							map(P->ctx.getGenotype(P.getId())).
							filter(this.acceptFilteredGenotype).
							anyMatch(g->(!(g==null ||
									!g.isCalled() ||
									!g.isAvailable() ||
									g.isNoCall() ||
									g.isHomRef())))
									;
					
									
					if(!in_pedigree) {
						if(CtxWriterFactory.this.dicardVariant) return;
						final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						vcb.filter(this.filter.getID());
						super.add(vcb.make());
						}
					else
						{
						boolean is_singleton;
						if(CtxWriterFactory.this.singleton!=IGNORE_SINGLETON) {
							is_singleton = true;
							for(final Allele alt:ctx.getAlternateAlleles()) {
								if(this.individuals.stream().
										map(P->ctx.getGenotype(P.getId())).
										filter(this.acceptFilteredGenotype).
										filter(g->g!=null && g.isCalled()&& g.getAlleles().contains(alt)).
										count() > CtxWriterFactory.this.singleton) 
									{
									is_singleton =false;
									break;
									}
								}
							}
						else
							{
							is_singleton=false;
							}
						if(is_singleton) {
							if(CtxWriterFactory.this.dicardVariant) return;
							final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
							vcb.filter(this.singletonFilter.getID());
							super.add(vcb.make());
							}
						else
							{
							super.add(ctx);
							}
						}
					}
				}
						
			@Override
			public VariantContextWriter open(final VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			
			}
	
	
	public VcfFilterNotInPedigree()
		{
		}
	 
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator in,
			final VariantContextWriter delegate
			) {
		final VariantContextWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progess = new SAMSequenceDictionaryProgress(in.getHeader()).logger(LOG);
		out.writeHeader(in.getHeader());
		while(in.hasNext())
			{
			out.add(progess.watch(in.next()));
			}
		progess.finish();
		out.close();
		return 0;
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, outputFile);
			}
		finally
			{
			CloserUtil.close(this.component);
			}
	}
	
	
	public static void main(final String[] args)
		{
		new VcfFilterNotInPedigree().instanceMainWithExit(args);
		}
	}
