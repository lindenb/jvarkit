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

*/
package com.github.lindenb.jvarkit.tools.gtf;

import java.io.BufferedReader; 
import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.iterator.LineIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3FeatureImpl;

/**

BEGIN_DOC

```
```

when the GTF contains a gene, all sub-section MUST be re-mapped.

```

```

END_DOC
*/
@Program(
		name="gffliftover",
		description="LiftOver GFF file.",
		creationDate="20201201",
		modificationDate="20201201",
		keywords= {"gff","liftover"}
		)
public class GffLiftOver
	extends Launcher
	{
	private static final Logger LOG = Logger.build(GffLiftOver.class).make();
	private LiftOver liftOver = null;
	
	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-f","--chain"},description="LiftOver file.",required=true)
	private File liftOverFile = null;
	@Parameter(names={"-x","--failed"},description="write  failing the liftOver here. Optional.")
	private Path failedFile = null;
	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;
	@Parameter(names={"-H","--header"},description= "include gtf header")
	private boolean include_header = false;

	private Interval applyLiftOver(final Gff3Feature feature) {
			return this.liftOver.liftOver(new Interval(
					feature.getContig(),
					feature.getStart(),
					feature.getEnd(),
					feature.getStrand().equals(Strand.NEGATIVE),
					feature.getID()
					),this.userMinMatch);
		}

	private boolean canLiftOver(final Gff3Feature feature) {
		return applyLiftOver(feature)!=null;
		}
	
	@Override
	public int doWork(final List<String> args)  {
		PrintWriter fail = null;
		try {
			fail = this.failedFile==null?
					new PrintWriter(new NullOuputStream()):
					super.openPathOrStdoutAsPrintWriter(this.failedFile)
					;			
			this.liftOver = new LiftOver(this.liftOverFile);			
			
			
			try(BufferedReader in = openBufferedReader(oneFileOrNull(args))) {
				final Gff3Codec codec  = new Gff3Codec();
				final LineIterator liter = new LineIterator(in);
				final FeatureCodecHeader header = codec.readHeader(liter);
				try(Gff3Writer p = new Gff3Writer(super.openPathOrStdoutAsStream(this.outputFile))) {
					while(!codec.isDone(liter)) {
						final Gff3Feature feature = codec.decode(liter);
						feature.getChildren().stream().filter(F->F.getType().equals("gene")).forEach(GENE->{
							if(!canLiftOver(GENE)) {
								continue;
								}
							p.addFeature(GENE);
							GENE.getChildren().stream().forEach(TRANSCRIPT->{
								final List<Gff3Feature> features = new ArrayList<>();
								features.add(TRANSCRIPT);
								TRANSCRIPT.getChildren().forEach(FEAT->features.add(FEAT));
								if(!features.stream().allMatch(F2->canLiftOver(F2))) {
									continue;
									}
								for(Gff3Feature feat: features) {
									final Interval lifted = applyLiftOver(feat);
									final Gff3FeatureImpl feat2 = new Gff3FeatureImpl(
											lifted.getContig(),
											feat.getSource(),
											feat.getType(),
											lifted.getStart(),
											lifted.getEnd(),
											feat.getScore(),
											lifted.getStrand(),
											feat.getPhase(),
											feat.getAttributes()
											);
									p.addFeature(feat2);
									}
								});
							
						
							});
						}
					codec.close(liter);
					p.flush();
					}
				}

			fail.flush();
			fail.close();
			fail = null;
			
			return RETURN_OK;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(fail);
			CloserUtil.close(p);
			}
		}
	
	public static void main(final String[] args)
		{
		new GffLiftOver().instanceMainWithExit(args);
		}
	}
