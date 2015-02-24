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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.knime.KnimeApplication;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
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
public class VcfIn extends AbstractCommandLineProgram
	implements KnimeApplication
	{
	private File outputFile=null;
	private boolean inverse=false;
	private boolean databaseIsTabix=true;
	private boolean userAltInDatabase=false;
	private int count_variants=0;
	public VcfIn()
		{
		}
	
	public int getCountCariant() {
		return count_variants;
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
		return "https://github.com/lindenb/jvarkit/wiki/VcfIn";
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
	
	@Override
	public void setOutputFile(File out) {
		this.outputFile=out;
		}
	
	private boolean sameChromPosRef(VariantContext ctx1,VariantContext ctx2)
		{
		if(!ctx1.getChr().equals(ctx2.getChr())) return false;
		if(ctx1.getStart()!=ctx2.getStart())return false;
		if(!ctx1.getReference().equals(ctx2.getReference()))return false;
		return true;
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
		VariantContextWriter w;
		if(this.outputFile==null)
			{
			w = VCFUtils.createVariantContextWriterToStdout();
			}
		else
			{
			w = VCFUtils.createVariantContextWriter(this.outputFile);
			}
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
		VcfIterator databaseIn=null;
		VcfIterator userVcfIn=null;
		try
			{
			databaseIn = VCFUtils.createVcfIterator(databaseVcfUri);
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
			
			vcw = this.openWriter(header);
			SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(userVcfDict);
			
			List<VariantContext> userContext = new ArrayList<VariantContext>();
			List<VariantContext> dbContext = new ArrayList<VariantContext>();
			while(userVcfIn.hasNext())
				{
				//fill both contextes
				
				if(dbContext.isEmpty() && databaseIn.hasNext())
					{
					VariantContext ctx= databaseIn.next();
					//not in user dict anyway
					if( userVcfDict.getSequence(ctx.getChr())==null)
						{
						dbContext.clear();
						continue;
						}
					
					dbContext.add(ctx);
					while(databaseIn.hasNext())
						{
						ctx = databaseIn.peek();
						if(!sameChromPosRef(ctx, dbContext.get(0)))
							{
							break;
							}
						dbContext.add(databaseIn.next());
						}
					}
				
				if(userContext.isEmpty() && userVcfIn.hasNext())
					{
					userContext.add(progress.watch(userVcfIn.next()));
					while(userVcfIn.hasNext())
						{
						VariantContext ctx = userVcfIn.peek();
						if(!sameChromPosRef(ctx, userContext.get(0)))
							{
							break;
							}
						userContext.add(progress.watch(userVcfIn.next()));
						}
					}
				
				if(userContext.isEmpty()) break;
				if(dbContext.isEmpty())
					{
					if(isInverse())
						{
						//flush buffer
						for(VariantContext ctx:userContext)
							{
							vcw.add(ctx);
							++count_variants;
							}
						//flush stream
						while(userVcfIn.hasNext())
							{
							vcw.add(progress.watch(userVcfIn.next()));
							++count_variants;
							}
						
						}
					//we're done
					break;
					}
				
				int i = userVcfDict.getSequence(  dbContext.get(0).getChr()).getSequenceIndex() -
						userVcfDict.getSequence(userContext.get(0).getChr()).getSequenceIndex()
						;
				if(i==0)
					{
					i = dbContext.get(0).getStart() - userContext.get(0).getStart();
					}
				if(i==0)
					{
					i = dbContext.get(0).getReference().compareTo(userContext.get(0).getReference());
					}
				//database before user
				if(i<0)
					{
					dbContext.clear();
					continue;
					}
				//database AFTER user
				if(i>0)
					{
					if(isInverse())
						{
						for(VariantContext ctx:userContext)
							{
							vcw.add(ctx);
							++count_variants;
							}
						}
					userContext.clear();
					continue;
					}
				
				for(VariantContext userCtx :userContext)
					{
					boolean keep=false;
					for(VariantContext dbCtx :dbContext)
						{
						if(!allUserAltFoundInDatabase(userCtx, dbCtx)) continue;
						keep=true;
						break;
						}
					if(isInverse()) keep=!keep;
					if(keep)
						{
						vcw.add(userCtx);
						++count_variants;
						}
					}
				userContext.clear();
				dbContext.clear();
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
			CloserUtil.close(databaseIn);
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
			this.count_variants=0;
			info("opening "+databaseVcfUri+" as tabix");
			tabix =  new TabixVcfFileReader(databaseVcfUri);
			in2 = (userVcfUri==null?VCFUtils.createVcfIteratorStdin():
				VCFUtils.createVcfIterator(userVcfUri));
			VCFHeader header1= in2.getHeader();
			vcw = this.openWriter(header1);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1.getSequenceDictionary());
			
			while(in2.hasNext())
				{
				final VariantContext userCtx= progress.watch(in2.next());
				checkKnimeCancelled();
				Iterator<VariantContext> iter= tabix.iterator(userCtx.getChr(),
						Math.max(1,userCtx.getStart()-1),
						userCtx.getEnd()+1);
				boolean keep=false;
				while(iter.hasNext())
					{
					VariantContext dbctx= iter.next();
					if(!sameChromPosRef(userCtx,dbctx)) continue;
					if(!allUserAltFoundInDatabase(userCtx, dbctx)) continue;
					keep=true;
					break;
					}
				if(this.isInverse()) keep=!keep;
				if(!keep) continue;
				vcw.add(userCtx);
				++count_variants;
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
	public void checkKnimeCancelled() {
		
		}
	
	@Override
	public int initializeKnime() {
		return 0;
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
			error("Illegal number of arguments.");
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
	public void disposeKnime() {
		
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"itA"))!=-1)
			{
			switch(c)
				{
				case 'i': this.setInverse(true);break;
				case 't': this.setDatabaseIsTabix(true);break;
				case 'A': this.setUserAltInDatabase(true);break;
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
		List<String> L=new ArrayList<String>();
		for(int i=opt.getOptInd(); i< args.length;++i)
			L.add(args[i]);
		return executeKnime(L);
		}
	public static void main(String[] args) {
		new VcfIn().instanceMainWithExit(args);
	}
	}
