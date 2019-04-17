package com.github.lindenb.jvarkit.util.jcommander;

import com.github.lindenb.jvarkit.io.IOUtilsTest;
import com.github.lindenb.jvarkit.lang.SmartComparatorTest;
import com.github.lindenb.jvarkit.lang.StringUtilsTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.util.CounterTest;
import com.github.lindenb.jvarkit.util.bio.IntervalParserTest;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodecTest;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactoryTest;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIteratorTest;
import com.github.lindenb.jvarkit.util.iterator.FilterIteratorTest;
import com.github.lindenb.jvarkit.util.iterator.LineIteratorTest;
import com.github.lindenb.jvarkit.util.iterator.MergingIteratorTest;
import com.github.lindenb.jvarkit.util.log.ProgressFactoryTest;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparatorTest;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;

@AlsoTest({IOUtilsTest.class,VCFUtilsTest.class,StringUtilsTest.class,IntervalParserTest.class,CounterTest.class,BedLineCodecTest.class,
	ProgressFactoryTest.class,ContigDictComparatorTest.class,SmartComparatorTest.class,
	EqualRangeIteratorTest.class,
	FilterIteratorTest.class,
	LineIteratorTest.class,
	MergingIteratorTest.class,
	SamRecordFilterFactoryTest.class
	})
public class LauncherTest {

}
