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
package com.github.lindenb.jvarkit.tools.bigwigtview;

import java.io.PrintStream;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.ansi.AnsiUtils;
import com.github.lindenb.jvarkit.ansi.AnsiUtils.AnsiColor;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.wig.BigWigReader;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
/**
BEGIN_DOC

## Motivation

plot set of bigwig files in terminal

## EXAMPLES

## Input

intut is one or more bigwigfile or a file with the '.list' suffix containing the path to the bigwigs
 
### Example

```
$ java -jar dist/jvarkit.jar bigwigtview --plain -region "chr1:1014238-1014330" src/test/resources/Uniqueness35bp.bigWig 


> src/test/resources/Uniqueness35bp.bigWig chr1:1014238-1014330
  1.050000 |                                                                                                     
  0.997500 |                                                               ######### #########             ##### 
  0.945000 |                                                               ######### #########             ##### 
  0.892500 |                                                               ######### #########             ##### 
  0.840000 |                                                               ######### #########             ##### 
  0.787500 |                                                               ######### #########             ##### 
  0.735000 |                                                               ######### #########             ##### 
  0.682500 |                                                               ######### #########             ##### 
  0.630000 |                                                               ######### #########             ##### 
  0.577500 |                                                               ######### #########             ##### 
  0.525000 |                                                               ######### #########             ##### 
  0.472500 |                                ###                            ######### ############# ############# 
  0.420000 |                                ###                            ######### ############# ############# 
  0.367500 |                                ###                            ######### ############# ############# 
  0.315000 |                                ###              ######### ############# ############# ############# 
  0.262500 |                                ###              ######### ############# ############# ############# 
  0.210000 |  ##                            ########### ############## ############# ############# ############# 
  0.157500 |  ##                            ########### ############## ############# ############# ############# 
  0.105000 |  ##                            ########### ############## ############# ############# ############# 
  0.052500 |  ##                            ########### ############## ############# ############# ############# 
  0.000000 | ####################################################################################################
      chr1 | 1014238 1014245 1014252 1014260 1014267 1014275 1014282 1014290 1014297 1014304 1014312 1014319 
```


END_DOC
*/
@Program(
		name="bigwigtview",
		description="view bigwig file coverage in a terminal",
		keywords={"wig","bigwig"},
		creationDate="20240704",
		modificationDate="20240704",
		jvarkit_amalgamion =  true
		)
public class BigWigTView extends Launcher {
	private static final Logger LOG = Logger.of(BigWigTView.class);
	
	@Parameter(names = { "-o", "--out" }, description = "Output is a setfile. "+ OPT_OUPUT_FILE_OR_STDOUT)
	private Path output = null;
	@Parameter(names = { "--interval","--region"}, description = "interval chr:start-end.",required = true)
	private String interval="";
	@Parameter(names = { "--height"}, description = "screen height")
	private int height=20;
	@Parameter(names = { "--width"}, description = "screen width")
	private int width=100;
	@Parameter(names = { "--plain"}, description = "disable ansi")
	private boolean disable_ansi=false;
	@Parameter(names = { "--max"}, description = "when binning data, use max value instead of average")
	private boolean use_max_value = false;
	@Override
	public int doWork(final List<String> args) {
		try {
			final List<String> wigs = IOUtils.unrollStrings(args);
			if(wigs.isEmpty()) {
				LOG.error("No bigwig defined");
				return -1;
				}
			if(width<10) {
				LOG.warning("width is too low. adjusting...");
				width=10;
				}
			if(height<5) {
				LOG.warning("height is too low. adjusting...");
				height=5;
				}
			
			final Locatable location = new IntervalParser().
					apply(this.interval).
					orElse(null);
			if(location==null) {
				LOG.error("cannot parse location from "+this.interval);
				return -1;
				}
			final double[] cov =new double[width];
			final int[] count =new int[width];
			try(PrintStream out=super.openPathOrStdoutAsPrintStream(output)) {
				Arrays.fill(cov, 0.0);
				Arrays.fill(count, 0);
				for(final String wig: wigs) {
					try( BigWigReader reader=new BigWigReader(wig)) {
						try(CloseableIterator<BigWigReader.WigItem> iter= reader.query(location)) {
							while(iter.hasNext()) {
								final BigWigReader.WigItem item=iter.next();
								for(int x=item.getStart();x<=item.getEnd();++x) {
									final int idx = (int)(((x-location.getStart())/(double)location.getLengthOnReference())*this.width);
									if(idx<0 || idx>=cov.length) continue;
									if(use_max_value) {
										cov[idx]=Math.max(item.getValue(),cov[idx]);
										}
									else
										{
										cov[idx]+=item.getValue();
										count[idx]++;
										}
									}
								}
							}
						//get average under point
						for(int i=0;i< cov.length && !this.use_max_value ;++i) {
							if(count[i]==0) continue;
							cov[i]/=count[i];
							}
						double maxValue = Arrays.stream(cov).max().orElse(1.0);
						if(maxValue==0.0) maxValue=1.0;
						maxValue+=maxValue*0.05;
						if(!disable_ansi) out.print(AnsiColor.CYAN.begin());
						out.println("> "+wig+ " "+interval);
						if(!disable_ansi) out.print(AnsiColor.CYAN.end());
						for(int y=this.height;y>=0;y--) {
							final double v= (y/(double)this.height)*maxValue;
							if(!disable_ansi) out.print(AnsiColor.GREEN.begin());
							out.print(String.format("%10f",v));
							if(!disable_ansi) out.print(AnsiColor.GREEN.end());
							out.print(" | ");
							if(!disable_ansi) out.print(AnsiColor.YELLOW.begin());
							for(int x=0;x< this.width;++x) {
								if(cov[x]<v) 
									{
									out.print(" ");
									}
								else if(!disable_ansi)
									{
									out.print(AnsiUtils.getHistogram(1f));
									}
								else
									{
									out.print("#");
									}
								}
							if(!disable_ansi) out.print(AnsiColor.YELLOW.end());
							out.println();
							}
						String txt = String.format("%10s",location.getContig());
						if(!disable_ansi) txt = AnsiColor.GREEN.colorize(txt);
						out.print(txt);
						out.print(" | ");
						int x=0;
						while(x< this.width) {
							int loc = location.getStart()+ (int)((x/(double)this.width)*location.getLengthOnReference());
							final String locstr= String.valueOf(loc);
							if(x+locstr.length()>this.width) break;
							if(!disable_ansi) out.print(AnsiColor.GREEN.begin());
							out.print(locstr);
							if(!disable_ansi) out.print(AnsiColor.GREEN.end());
							out.print(" ");
							x+=locstr.length()+1;
							}
						out.println();
						out.println();
						}
					}
				out.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new BigWigTView().instanceMainWithExit(args);
	}

}
