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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;

/**
 * VcfIn
 *
 */
public class VcfIn extends AbstractVCFFilter3
	{
	private boolean inverse=false;
	private boolean databaseIsTabix=false;
	private boolean userAltInDatabase=false;
	public VcfIn()
		{
		}
	
	/**
	 * print user variant only if not found in database VCF
	 * @param b
	 */
	public void setInverse(boolean b)
		{
		this.inverse=b;
		}
	
	/** print user variant that are not part of database ? */
	public boolean isInverse()
		{
		return inverse;
		}
	
	public void setDatabaseIsTabix(boolean databaseIsTabix) {
		this.databaseIsTabix = databaseIsTabix;
		}
	
	public void setUserAltInDatabase(boolean userAltInDatabase) {
		this.userAltInDatabase = userAltInDatabase;
		}
	
	public boolean isUserAltInDatabase() {
		return userAltInDatabase;
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VcfIn";
		}

		@Override
	public String getProgramDescription() {
		return "Only prints variants that are contained/not contained into another VCF.";
		}
	
		
		
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -i : inverse. Print variant that are not part of the VCF-database.");
		out.println(" -t : database is tabix-ed.");
		out.println(" -A : ALL user ALT must be found in VCF-database ALT.");
		super.printOptions(out);
		}
	
	
	
	private boolean allUserAltFoundInDatabase(
			VariantContext userVariants,
			VariantContext databaseVariants
			)
		{
		if(!isUserAltInDatabase()) return true;
		Set<Allele> user_alts=new HashSet<Allele>(userVariants.getAlternateAlleles());
		user_alts.removeAll(databaseVariants.getAlternateAlleles());
		return user_alts.isEmpty();
		}
	
	private VariantContextWriter openWriter(VCFHeader header)
		throws IOException
		{
		VariantContextWriter w = super.createVariantContextWriter();
		VCFHeader h2= new VCFHeader(header);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		w.writeHeader(h2);
		return w;
		}
	
	private int scanFileSorted(String databaseVcfUri,String userVcfUri)
		{
		VariantContextWriter vcw=null;
		VcfIterator userVcfIn=null;
		EqualRangeVcfIterator equalRangeDbIter=null;
		try
			{
			userVcfIn = (userVcfUri==null?VCFUtils.createVcfIteratorStdin():
				VCFUtils.createVcfIterator(userVcfUri));
			VCFHeader header = userVcfIn.getHeader();
			SAMSequenceDictionary userVcfDict = new VCFHeader(header).getSequenceDictionary();
			/// NO need if(dict1==null)
			if(userVcfDict==null)
				{
				error("NO SAM sequence Dict in user VCF "+(userVcfUri==null?"stdin":userVcfUri));
				return -1;
				}
			Comparator<VariantContext> vcfComparator =
					VCFUtils.createTidPosComparator(userVcfDict)
					;
			equalRangeDbIter = new EqualRangeVcfIterator(
					VCFUtils.createVcfIterator(databaseVcfUri),vcfComparator);

			vcw = this.openWriter(header);
			SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(userVcfDict);
			
			
			
			while(userVcfIn.hasNext())
				{
				VariantContext ctx = progress.watch(userVcfIn.next());
				//fill both contextes
				List<VariantContext> dbContexes = new ArrayList<VariantContext>(equalRangeDbIter.next(ctx));
				
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
				
				boolean keep=!dbContexes.isEmpty();
				
				if(isInverse()) keep=!keep;
				if(keep)
					{
					vcw.add(ctx);
					incrVariantCount();
					}
					
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(equalRangeDbIter);
			CloserUtil.close(userVcfIn);
			CloserUtil.close(vcw);
			}
		}
	
	private int scanUsingTabix(String databaseVcfUri,String userVcfUri)
		{
		VariantContextWriter vcw=null;
		TabixVcfFileReader tabix=null;
		VcfIterator in2=null;
		try
			{
			setVariantCount(0);
			info("opening "+databaseVcfUri+" as tabix");
			tabix =  new TabixVcfFileReader(databaseVcfUri);
			in2 = (userVcfUri==null?VCFUtils.createVcfIteratorStdin():
				VCFUtils.createVcfIterator(userVcfUri));
			VCFHeader header1= in2.getHeader();
			vcw = this.openWriter(header1);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1.getSequenceDictionary());
			
			while(in2.hasNext() && !this.checkOutputError())
				{
				final VariantContext userCtx= progress.watch(in2.next());
				checkKnimeCancelled();
				Iterator<VariantContext> iter= tabix.iterator(userCtx.getContig(),
						Math.max(1,userCtx.getStart()-1),
						userCtx.getEnd()+1);
				boolean keep=false;
				while(iter.hasNext())
					{
					VariantContext dbctx= iter.next();
					if(! dbctx.getContig().equals(userCtx.getContig())) continue;
					if(dbctx.getStart()!=userCtx.getStart()) continue;
					if(! dbctx.getReference().equals(userCtx.getReference())) continue;
					if(!allUserAltFoundInDatabase(userCtx, dbctx)) continue;
					keep=true;
					break;
					}
				if(this.isInverse()) keep=!keep;
				if(!keep) continue;
				vcw.add(userCtx);
				incrVariantCount();
				}
			progress.finish();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(tabix);
			CloserUtil.close(in2);
			CloserUtil.close(vcw);
			}
		}
	
	
	@Override
	public int executeKnime(List<String> args)
		{
		String databaseVcfUri = null;
		String userVcfUri = null;

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
			error(getMessageBundle("illegal.number.of.arguments"));
			return -1;
			}

		if(this.databaseIsTabix)
			{
			return this.scanUsingTabix(databaseVcfUri, userVcfUri);
			}
		else
			{	
			return scanFileSorted(databaseVcfUri, userVcfUri);
			}
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"itAo:"))!=-1)
			{
			switch(c)
				{
				case 'i': this.setInverse(true);break;
				case 't': this.setDatabaseIsTabix(true);break;
				case 'A': this.setUserAltInDatabase(true);break;
				case 'o': this.setOutputFile(opt.getOptArg()); break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}
	public static void main(String[] args) {
		new VcfIn().instanceMainWithExit(args);
	}
	}
