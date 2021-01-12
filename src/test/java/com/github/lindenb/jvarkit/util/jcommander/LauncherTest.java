package com.github.lindenb.jvarkit.util.jcommander;

import com.github.lindenb.jvarkit.io.IOUtilsTest;
import com.github.lindenb.jvarkit.lang.SmartComparatorTest;
import com.github.lindenb.jvarkit.lang.StringUtilsTest;
import com.github.lindenb.jvarkit.lang.primitive.DoubleArrayTest;
import com.github.lindenb.jvarkit.lang.primitive.IntArrayTest;
import com.github.lindenb.jvarkit.samtools.reference.TwoBitSequenceFileTest;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactoryTest;
import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.math.DiscreteMedianTest;
import com.github.lindenb.jvarkit.util.CounterTest;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodecTest;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactoryTest;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIteratorTest;
import com.github.lindenb.jvarkit.util.iterator.FilterIteratorTest;
import com.github.lindenb.jvarkit.util.iterator.LineIteratorTest;
import com.github.lindenb.jvarkit.util.iterator.MergingIteratorTest;
import com.github.lindenb.jvarkit.util.log.ProgressFactoryTest;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparatorTest;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;

@AlsoTest({IOUtilsTest.class,VCFUtilsTest.class,StringUtilsTest.class,IntervalParserFactoryTest.class,CounterTest.class,BedLineCodecTest.class,
	ProgressFactoryTest.class,ContigDictComparatorTest.class,SmartComparatorTest.class,
	EqualRangeIteratorTest.class,
	FilterIteratorTest.class,
	LineIteratorTest.class,
	MergingIteratorTest.class,
	SamRecordFilterFactoryTest.class,
	TwoBitSequenceFileTest.class,
	DiscreteMedianTest.class,
	IntArrayTest.class,
	DoubleArrayTest.class
	})
public class LauncherTest {

}
