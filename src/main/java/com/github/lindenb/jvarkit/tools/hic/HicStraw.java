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
package com.github.lindenb.jvarkit.tools.hic;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.hic.HicReader;
import com.github.lindenb.jvarkit.hic.HicReaderFactory;
import com.github.lindenb.jvarkit.hic.Normalization;
import com.github.lindenb.jvarkit.hic.Unit;
import com.github.lindenb.jvarkit.io.CustomSeekableStreamFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;

/**
BEGIN_DOC

### Example


```

END_DOC
 */

@Program(name="hicstraw",
	description="Query a Hi-C file",
	keywords={"hic"},
	creationDate="20190613",
	modificationDate="20190613",
	generate_doc=false
	)
public class HicStraw  extends Launcher {
	private static final Logger LOG = Logger.build(HicFileInfo.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-i","--interval1"},description="Interval 1",required=true)
	private String interval1Str = null;
	@Parameter(names={"-j","--interval2"},description="Interval 2",required=true)
	private String interval2Str = null;
	@Parameter(names={"-u","--unit"},description="Unit")
	private Unit unit = Unit.BP;
	@Parameter(names={"-n","--normalization"},description="normalization")
	private Normalization norm = Normalization.NONE;
	@Parameter(names={"-b","--bin"},description="bin size")
	private int binSize = 1;

	private final HicReader.QueryCallBack defaultCallBack = new HicReader.QueryCallBack() {
		};
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			final ISeekableStreamFactory seekableStreamFactory = new CustomSeekableStreamFactory();
			
			for(final String input :args) {
			
				
				try(final HicReader hicReader = new HicReaderFactory().
							setSeekableStreamFactory(seekableStreamFactory).
							open(input)) { 
				
					final Function<String,Locatable > parseInterval = (S)->{
						final Optional<Locatable> loc = hicReader.parseInterval(S);
						if(!loc.isPresent()) {
							LOG.error("bad interval : \""+S+"\" available are "+ hicReader.
									getDictionary().getSequences().stream().
									map(SR->SR.getSequenceName()).collect(Collectors.joining(" ; ")));
							return null;
							}
						return loc.get();
						};
					
					if(!hicReader.getBasePairResolutions().contains(this.binSize)) {
						LOG.error("bad binSize : \""+this.binSize+"\" available are "+ hicReader.getBasePairResolutions().stream().map(S->String.valueOf(S)).collect(Collectors.joining(" ; ")));
						return -1;
						}
						
					final Locatable loc1 = parseInterval.apply(this.interval1Str);
					if(loc1==null) return -1;
					final Locatable loc2 = parseInterval.apply(this.interval2Str);
					if(loc2==null) return -1 ;
					hicReader.query(loc1, loc2,norm, this.binSize, this.unit, this.defaultCallBack);
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args) {
		new HicStraw().instanceMainWithExit(args);
		}
}
