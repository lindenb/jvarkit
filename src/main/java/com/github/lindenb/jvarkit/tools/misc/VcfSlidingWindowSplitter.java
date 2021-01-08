/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.ArrayList;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

/**

BEGIN_DOC

### Example

```
$ java -jar dist/vcfwindowsplitter.jar -n 2 -w 1000 -s 500 -o jeter.zip -m jeter.manifest src/test/resources/rotavirus_rf.vcf.gz 
[INFO][VcfSlidingWindowSplitter]. Completed. N=45. That took:0 second

$ head jeter.manifest  | column -t
#chrom  start  end   window          path                                                     Count_Variants
RF02    250    877   RF02:0-1000     35/a02abf73216a95f62458f41dfbbf79/RF02_0_1000.vcf.gz     3
RF02    577    877   RF02:500-1500   b5/ed3c205185d506e382c437cc010dee/RF02_500_1500.vcf.gz   2
RF02    1725   1965  RF02:1000-2000  41/026109d82e6f553dd0b6244ceffe5a/RF02_1000_2000.vcf.gz  2
RF02    1725   1965  RF02:1500-2500  2f/82d6416c74c8b2427fb71a696bf491/RF02_1500_2500.vcf.gz  2
RF03    1220   1242  RF03:500-1500   3b/11b08845e1d70e4b4c6ab3a5a279ee/RF03_500_1500.vcf.gz   2
RF03    1220   1708  RF03:1000-2000  dc/e6069ebdb06ef2cfd85e96331d29d2/RF03_1000_2000.vcf.gz  4
RF03    1687   2315  RF03:1500-2500  af/b0daefa93f1423494ed6da29ab49c7/RF03_1500_2500.vcf.gz  5
RF03    2149   2573  RF03:2000-3000  f7/ffe3fb2204dd071e7d40d8e0e1151e/RF03_2000_3000.vcf.gz  4
RF04    886    991   RF04:0-1000     5b/9e9699f08473c9fa7cd550247c36ff/RF04_0_1000.vcf.gz     2

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     1359  2019-06-19 14:48   35/a02abf73216a95f62458f41dfbbf79/RF02_0_1000.vcf.gz
     1283  2019-06-19 14:48   b5/ed3c205185d506e382c437cc010dee/RF02_500_1500.vcf.gz
     1301  2019-06-19 14:48   41/026109d82e6f553dd0b6244ceffe5a/RF02_1000_2000.vcf.gz
     1301  2019-06-19 14:48   2f/82d6416c74c8b2427fb71a696bf491/RF02_1500_2500.vcf.gz
     1281  2019-06-19 14:48   3b/11b08845e1d70e4b4c6ab3a5a279ee/RF03_500_1500.vcf.gz
     1392  2019-06-19 14:48   dc/e6069ebdb06ef2cfd85e96331d29d2/RF03_1000_2000.vcf.gz
     1462  2019-06-19 14:48   af/b0daefa93f1423494ed6da29ab49c7/RF03_1500_2500.vcf.gz
     1418  2019-06-19 14:48   f7/ffe3fb2204dd071e7d40d8e0e1151e/RF03_2000_3000.vcf.gz
     1265  2019-06-19 14:48   5b/9e9699f08473c9fa7cd550247c36ff/RF04_0_1000.vcf.gz
     1405  2019-06-19 14:48   4c/75302cbe2cd7f16759c5e864031e03/RF04_500_1500.vcf.gz
     1482  2019-06-19 14:48   d4/7c942d8e574e6bdfccafbeb767675e/RF04_1000_2000.vcf.gz
     1356  2019-06-19 14:48   7f/af3063f079076117e4186fab17768b/RF04_1500_2500.vcf.gz
     1414  2019-06-19 14:48   e6/73552f774e54d0c0bd8bfc31192335/RF05_0_1000.vcf.gz
     1410  2019-06-19 14:48   95/3b7f4b269e5dc25550ed7610c4be52/RF05_500_1500.vcf.gz
     1274  2019-06-19 14:48   1e/a8b7283d053ccf22f44301504d522f/RF05_1000_2000.vcf.gz
     1427  2019-06-19 14:48   ea/ef4a24f964402115bfbf6b9f25253e/RF06_0_1000.vcf.gz
     1504  2019-06-19 14:48   d5/487160765c2bc2a1d286e534c41d06/RF06_500_1500.vcf.gz
     1406  2019-06-19 14:48   15/c37e21e6933c610a8bf53d833fcea9/RF07_0_1000.vcf.gz
     1270  2019-06-19 14:48   ff/af5bab3250e5ea8106efc4ca424e4f/RF07_500_1500.vcf.gz
     1269  2019-06-19 14:48   4b/9c386cb83786e3a61bbc54b949d5f9/RF08_0_1000.vcf.gz
     1270  2019-06-19 14:48   bb/fc42ed842a1d5f76bd19986655c7f7/RF08_500_1500.vcf.gz
     1333  2019-06-19 14:48   07/9587be0df182f9e938d412ffd9ffa6/RF09_0_1000.vcf.gz
     1341  2019-06-19 14:48   b1/5a51dc37c5b128a2edd36a68fcb9cb/RF10_0_1000.vcf.gz
---------                     -------
    31223                     23 files

```

END_DOC
*/
@Program(
		name="vcfwindowsplitter",
		description="Split VCF by sliding window",
		creationDate="20190619",
		modificationDate="20191129",
		keywords= {"vcf","sliding","window"}
		)
public class VcfSlidingWindowSplitter
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfSlidingWindowSplitter.class).make();
		
	@Parameter(names={"-o","--output"},description= ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-m","--manifest"},description="Manifest Bed file output containing chrom/start/end of each gene")
	private Path manifestFile = null;
	@Parameter(names={"-w","-W","--window-size"},description="Sliding window size. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int window_size = 1_000_000;
	@Parameter(names={"-s","-S","--window-shift"},description="Sliding window shift. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int window_shift = 500_000;
	@Parameter(names={"-n","--min-variant"},description="Minimum number of variants required to write a vcf. don't write if num(variant) < 'x' ")
	private int min_number_of_ctx = 1;
	@Parameter(names={"-M","--max-variant"},description="Maximum number of variants required to write a vcf. don't write if num(variant) > 'x' . '<=0' is ignore")
	private int max_number_of_ctx = -1;


	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	private static class WinAndLine {
		final Interval interval;
		final String ctx;
		WinAndLine(final Interval interval,final String ctx) {
			this.interval = interval;
			this.ctx = ctx;
		}
	}
	
	private static class WinAndLineComparator1
		implements Comparator<WinAndLine>
		{
		@Override
		public int compare(final WinAndLine o1, final WinAndLine o2) {
			return o1.interval.compareTo(o2.interval);
			}
		}
	private static class WinAndLineComparator2
		extends WinAndLineComparator1
		{
		final CharSplitter tab = CharSplitter.TAB;
		@Override
		public int compare(final WinAndLine o1, final WinAndLine o2) {
			int i = super.compare(o1, o2);
			if(i!=0) return i;
			final String tokens1[] = this.tab.split(o1.ctx,5);
			final String tokens2[] = this.tab.split(o2.ctx,5);
			if(!tokens1[0].equals(tokens2[0])) {
				throw new IllegalStateException("not same contig???");
			}
			i = Integer.parseInt(tokens1[1]) - Integer.parseInt(tokens2[1]);
			if(i!=0) return i;
			i = Allele.create(tokens1[3],true).compareTo( Allele.create(tokens2[3],true));
			return i;
			}
		}
	
	private static class WinAndLineCodec extends AbstractDataCodec<WinAndLine>
		{
		@Override
		public WinAndLine decode(final DataInputStream dis) throws IOException {
			String ctg;
			try {
				ctg = dis.readUTF();
			} catch(final EOFException err) { return null;}
			final int start = dis.readInt();
			final int end = dis.readInt();
			final String line = AbstractDataCodec.readString(dis);
			return new WinAndLine(new Interval(ctg,start,end),line);
		}
		@Override
		public void encode(final DataOutputStream dos,final  WinAndLine object) throws IOException {
			dos.writeUTF(object.interval.getContig());
			dos.writeInt(object.interval.getStart());
			dos.writeInt(object.interval.getEnd());
			AbstractDataCodec.writeString(dos, object.ctx);
			}
		@Override
		public AbstractDataCodec<WinAndLine> clone() {
			return new WinAndLineCodec();
			}
		}

		
	public VcfSlidingWindowSplitter()
		{
		}
	
	
	private int run(final List<String> args) {
		SortingCollection<WinAndLine> sortingcollection=null;
		BufferedReader in = null;
		FileOutputStream fos = null;
		CloseableIterator<WinAndLine> iter=null;
		ArchiveFactory archiveFactory = null;
		PrintWriter manifest = null;
		
		final Function<VariantContext,List<Interval>> makeWindows = (ctx)-> {
			final int ctx_start = ctx.getStart();
			final int ctx_end = ctx.getEnd();
			final List<Interval> list = new ArrayList<>();
		    int right = ctx_start -  ctx_start%this.window_shift;
		    while (right + this.window_size >= ctx_start )
		    	{
		    	right -= this.window_shift;
		    	}
		    
			while(right <= ctx_end)
		    		{
		    		final int left = right + this.window_size;

		    		if(  right>0 && CoordMath.overlaps(right,left,ctx_start,ctx_end) )  {
						list.add(new Interval(ctx.getContig(),right,left));
						}
		    		right += this.window_shift;
		    		}
		    	return list;
		    };
		
		try {
			final Path tmpVcf = Files.createTempFile("tmp.", ".vcf.gz");
			
			archiveFactory = ArchiveFactory.open(this.outputFile);
			
			manifest = new PrintWriter(this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(manifestFile));
			manifest.println("#chrom\tstart\tend\twindow\tpath\tCount_Variants");

			in = super.openBufferedReader(oneFileOrNull(args));
			final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(in);
			
			// read variants
			final ProgressFactory.Watcher<VariantContext> progess= ProgressFactory.newInstance().dictionary(cah.header).logger(LOG).build();
			String prevCtg = null;
			for(;;)
				{
				final String line = in.readLine();
				
				final VariantContext ctx = (line==null?null:cah.codec.decode(line));
				if(ctx!=null) progess.apply(ctx);
				
				if(ctx==null || !ctx.getContig().equals(prevCtg))
					{
					if(sortingcollection!=null) {
						sortingcollection.doneAdding();
						iter = sortingcollection.iterator();
						final EqualRangeIterator<WinAndLine> eqiter = new EqualRangeIterator<>(iter,
							new WinAndLineComparator1()
							);
						
						while(eqiter.hasNext()) {
							final List<WinAndLine> buffer = eqiter.next();
							if ( buffer.size() < this.min_number_of_ctx ) continue;
							if ( this.max_number_of_ctx > 0 && buffer.size() > this.max_number_of_ctx ) continue;
							final WinAndLine first = buffer.get(0);
							
							final VariantContextWriter out = VCFUtils.createVariantContextWriterToPath(tmpVcf);
							final VCFHeader header2=addMetaData(new VCFHeader(cah.header));
							header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+".interval",first.interval.getContig()+":"+first.interval.getStart()+"-"+first.interval.getEnd()));
							out.writeHeader(header2);
							int minPos=Integer.MAX_VALUE;
							int maxPos=0;
							for(final WinAndLine kl:buffer) {
								final VariantContext ctx2 = cah.codec.decode(kl.ctx);
								minPos = Math.min(ctx2.getStart(), minPos);
								maxPos = Math.max(ctx2.getEnd(), maxPos);
								out.add(ctx2);
								}
							out.close();
							
							final String md5 = StringUtils.md5(first.interval.getContig()+":"+first.interval.getStart()+"-"+first.interval.getEnd());
							final String filename =  md5.substring(0,2) + 
										File.separatorChar + 
										md5.substring(2) + 
										File.separator + 
										first.interval.getContig()+"_"+first.interval.getStart()+"_"+first.interval.getEnd() + 
										".vcf.gz";
							
							
							try( final OutputStream os = archiveFactory.openOuputStream(filename)) {
								IOUtils.copyTo(tmpVcf, os);
								os.flush();
								}
							
							manifest.print(prevCtg);
							manifest.print('\t');
							manifest.print(minPos-1);
							manifest.print('\t');
							manifest.print(maxPos);
							manifest.print('\t');
							manifest.print(first.interval.getContig()+":"+first.interval.getStart()+"-"+first.interval.getEnd());
							manifest.print('\t');
							manifest.print((archiveFactory.isTarOrZipArchive()?"":this.outputFile.toString()+File.separator)+filename);
							manifest.print('\t');
							manifest.println(buffer.size());
							}
						
						eqiter.close();
						iter.close();
						
					
						}
					sortingcollection = null;
					if(ctx==null) break;
					
					prevCtg = ctx.getContig();
					}
				
			 
				
		    	for(final Interval win: makeWindows.apply(ctx))
		    		{		    		
		    		if(sortingcollection==null) {
						sortingcollection = SortingCollection.newInstance(
								WinAndLine.class,
								new WinAndLineCodec(),
								new WinAndLineComparator2(),
								this.writingSortingCollection.getMaxRecordsInRam(),
								this.writingSortingCollection.getTmpPaths()
								);
						sortingcollection.setDestructiveIteration(true);
						}
					
				sortingcollection.add( new WinAndLine( win, line ) );
		    		}

				}
			progess.close();
			manifest.flush();
			manifest.close();
			archiveFactory.close();
			Files.deleteIfExists(tmpVcf);
			return RETURN_OK;
			}
		catch(final Exception err) 
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			if(sortingcollection!=null) sortingcollection.cleanup();
			CloserUtil.close(in);
			CloserUtil.close(fos);
			CloserUtil.close(manifest);
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
    	if(this.window_size<=0) {
    		LOG.error("Bad window size.");
    		return -1;
    	}
    	if(this.window_shift<=0) {
    		LOG.error("Bad window shift.");
    		return -1;
    	}
    	
    	if(this.min_number_of_ctx<=0) {
    		LOG.error("Bad minimum number of variants");
    		return -1;
    	}
    	
		try
			{
			return run(args);
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	 	
	
	public static void main(final String[] args)
		{
		new VcfSlidingWindowSplitter().instanceMainWithExit(args);
		}
	}
