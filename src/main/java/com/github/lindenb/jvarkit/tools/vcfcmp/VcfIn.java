/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;

/**
 * VcfIn
 *
 */
public class VcfIn extends AbstractVcfIn
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfIn.class);

	public VcfIn()
		{
		}
		
	
	private boolean allUserAltFoundInDatabase(
			final VariantContext userVariants,
			final VariantContext databaseVariants
			)
		{
		if(!super.isUserAltInDatabase()) return true;
		final Set<Allele> user_alts=new HashSet<Allele>(userVariants.getAlternateAlleles());
		user_alts.removeAll(databaseVariants.getAlternateAlleles());
		return user_alts.isEmpty();
		}
	
	@Override
	protected VCFHeader addMetaData(VCFHeader header) {
		if(!super.filterIn.isEmpty()) {
			header.addMetaDataLine(new VCFFilterHeaderLine(super.filterIn,
					"Variant overlapping database."));
			}
		if(!super.filterOut.isEmpty()) {
			header.addMetaDataLine(new VCFFilterHeaderLine(super.filterOut,
					"Variant non overlapping database."));
			}
		return super.addMetaData(header);
		}
	
	private void addVariant(final VariantContextWriter w,final VariantContext ctx,boolean keep)
		{
		if(isInverse()) keep=!keep;
		if(!this.filterIn.isEmpty())
			{
			if(keep){
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.filter(this.filterIn);
				w.add(vcb.make());
				}
			else
				{
				w.add(ctx);
				}
			}
		else  if(!this.filterOut.isEmpty()) {
			if(keep){
				w.add(ctx);
				}
			else
				{
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.filter(this.filterOut);
				w.add(vcb.make());
				}
			}
		else
			{
			if(keep) {
				w.add(ctx);
				} else
				{
					/* don't print */
				}
			}
		}
	
	/* public for knime */
	public Collection<Throwable> scanFileSorted(
			final VariantContextWriter vcw,
			final String databaseVcfUri,
			final VcfIterator userVcfIn
			)
		{
		EqualRangeVcfIterator equalRangeDbIter=null;
		try
			{
			final VCFHeader header = new VCFHeader(userVcfIn.getHeader());
			final SAMSequenceDictionary userVcfDict = header.getSequenceDictionary();
			/// NO need if(dict1==null)
			if(userVcfDict==null)
				{
				return wrapException("NO SAM sequence Dict in user VCF");
				}
			final Comparator<VariantContext> vcfComparator =
					VCFUtils.createTidPosComparator(userVcfDict)
					;
			equalRangeDbIter = new EqualRangeVcfIterator(
					VCFUtils.createVcfIterator(databaseVcfUri),vcfComparator);

			this.addMetaData(header);
			vcw.writeHeader(header);
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(userVcfDict);
			
			while(userVcfIn.hasNext())
				{
				final VariantContext ctx = progress.watch(userVcfIn.next());
				//fill both contextes
				final List<VariantContext> dbContexes = new ArrayList<VariantContext>(equalRangeDbIter.next(ctx));
				
				int i=0;
				while(i< dbContexes.size())
					{
					if( dbContexes.get(i).getReference().equals(ctx.getReference()) &&
						allUserAltFoundInDatabase(ctx, dbContexes.get(i)))
						{
						++i;
						}
					else
						{
						dbContexes.remove(i);
						}
					}
				
				final boolean keep=!dbContexes.isEmpty();
				addVariant(vcw,ctx,keep);
				if(vcw.checkError()) break;
				}
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(equalRangeDbIter);
			CloserUtil.close(userVcfIn);
			CloserUtil.close(vcw);
			}
		}
	/* public for knime */
	public Collection<Throwable> scanUsingTabix(final VariantContextWriter vcw,final String databaseVcfUri,final VcfIterator in2)
		{
		TabixVcfFileReader tabix=null;
		try
			{
			LOG.info("opening "+databaseVcfUri+" as tabix");
			tabix =  new TabixVcfFileReader(databaseVcfUri);
			final VCFHeader header1= new VCFHeader(in2.getHeader());
			this.addMetaData(header1);
			vcw.writeHeader(header1);
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1.getSequenceDictionary());
			
			while(in2.hasNext() && !vcw.checkError())
				{
				final VariantContext userCtx= progress.watch(in2.next());
				
				final Iterator<VariantContext> iter= tabix.iterator(userCtx.getContig(),
						Math.max(1,userCtx.getStart()-1),
						userCtx.getEnd()+1);
				boolean keep=false;
				while(iter.hasNext())
					{
					final VariantContext dbctx= iter.next();
					if(! dbctx.getContig().equals(userCtx.getContig())) continue;
					if(dbctx.getStart()!=userCtx.getStart()) continue;
					if(! dbctx.getReference().equals(userCtx.getReference())) continue;
					if(!allUserAltFoundInDatabase(userCtx, dbctx)) continue;
					keep=true;
					break;
					}
				
				addVariant(vcw,userCtx,keep);
				if(vcw.checkError()) break;
				}
			progress.finish();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(tabix);
			CloserUtil.close(in2);
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		throw new IllegalStateException("shouldn't be called");
		}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		if(!super.filterIn.isEmpty() && !super.filterOut.isEmpty()) {
			return wrapException("Option -"+OPTION_FILTERIN+"  and  -"+OPTION_FILTEROUT+" both defined.");
		}
		if(super.inverse && (!super.filterIn.isEmpty() || !super.filterOut.isEmpty())) {
			return wrapException("Option -"+OPTION_INVERSE+" cannot be used when Option -"+OPTION_FILTERIN+" or  -"+OPTION_FILTEROUT+" is defined.");
		}
		
		final List<String> args=super.getInputFiles(); 
		String databaseVcfUri;
		String userVcfUri;
		if(args.size()==1)
			{
			databaseVcfUri = args.get(0);
			userVcfUri =null;
			}
		else if(args.size()==2)
			{
			databaseVcfUri = args.get(0);
			userVcfUri = args.get(1);
			}
		else
			{
			return wrapException(getMessageBundle("illegal.number.of.arguments"));
			}

		VariantContextWriter w=null;
		VcfIterator in=null;
		try {
			in = (userVcfUri==null?
					VCFUtils.createVcfIteratorFromInputStream(stdin()):
					VCFUtils.createVcfIterator(userVcfUri)
					);
			w= super.openVariantContextWriter();
			if(super.databaseIsTabix)
				{
				return this.scanUsingTabix(w,databaseVcfUri, in);
				}
			else
				{
				return this.scanFileSorted(w,databaseVcfUri, in);
				}
			} catch (Exception err) {
				return wrapException(err);
			} finally
			{
			CloserUtil.close(in);
			CloserUtil.close(w);
			}
		
		
		}
	
	
	
	public static void main(String[] args) {
		new VcfIn().instanceMainWithExit(args);
	}
	}
