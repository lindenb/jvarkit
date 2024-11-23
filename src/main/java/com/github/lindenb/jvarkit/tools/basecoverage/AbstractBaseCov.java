/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.basecoverage;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.OptionalDouble;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;

public abstract class AbstractBaseCov extends Launcher {
	protected static final String FORMAT_NORM_DEPH= "MD";
	
	@Parameter(names={"--mapq"},description=" min mapping quality.")
	private int mapping_quality=1;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@ParametersDelegate
	protected WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	
	protected int[] getCoverage(final SamReader sr,final Locatable queryInterval) {
		System.gc();
		final int[] coverage_int = new int[queryInterval.getLengthOnReference()];
		Arrays.fill(coverage_int,0);
		try(CloseableIterator<SAMRecord> it= sr.queryOverlapping(queryInterval.getContig(), queryInterval.getStart(), queryInterval.getEnd())) {
			while(it.hasNext()) {
				final SAMRecord rec =it.next();
				if(!SAMRecordDefaultFilter.accept(rec, this.mapping_quality)) continue;
				for(AlignmentBlock ab : rec.getAlignmentBlocks()) {
					for(int x=0;x < ab.getLength();++x) {
						final int array_index = (ab.getReferenceStart()+x) - queryInterval.getStart();
						if(array_index<0 ) continue;
						if(array_index>=coverage_int.length) break;
						coverage_int[array_index]++;
						}
					}
				}
			}
		return coverage_int;
		}
	
	protected OptionalDouble getMedian(final int[] coverage_int) {
		final DiscreteMedian<Integer> dm = new DiscreteMedian<>();
		for(int array_index=0;array_index< coverage_int.length;++array_index) {
			dm.add(coverage_int[array_index]);
			}
		return dm.getMedian();
		}
	
	protected float[] normalizeOnMedian(final int[] coverage_int,double median_cov) {
		float[] coverage_norm = new float[coverage_int.length];	
		for(int array_index=0;array_index< coverage_int.length;++array_index) {
			coverage_norm[array_index]=(float)(coverage_int[array_index]/median_cov);
			}
		return coverage_norm;
		}
	
	protected AbstractBaseCov() {
		}

	}
