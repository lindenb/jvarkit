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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.cmpbams;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
/*
BEGIN_DOC

## Example

Question was : "mapping with BWA produced a variation one year ago. We then mapped the same fastq with two different sets of parameters, but we cannot find the variant anymore. Has the mapping changed ?"

Extract the read names from the original BAM:
```
samtools view  file1.bam K01:2179-2179 |\
 cut -d '	' -f 1  | sort | uniq > names.txt
```
 
Use [[SamGrep]] to retieve the reads in the 3 bams:

```
java -jar dist/samgrep.jar -f names.txt file1.bam &gt; tmp1.sam
java -jar dist/samgrep.jar -f names.txt file2.bam &gt; tmp2.sam
java -jar dist/samgrep.jar -f names.txt file3.bam &gt; tmp3.sam
```

Run CmpBams

```
$ java -jar dist/cmpbams.jar -F -C tmp1.sam tmp2.sam tmp3.sam

#READ-Name	tmp1.sam tmp2.sam|tmp1.sam tmp3.sam|tmp2.sam tmp3.sam	tmp1.sam	tmp2.sam	tmp3.sam
HWI-1KL149:20:C1CU7ACXX:1:1101:17626:32431/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:1:1101:17626:32431/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1102:16831:71728/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1102:16831:71728/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1105:3309:27760/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1105:3309:27760/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1106:2914:12111/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:1:1106:2914:12111/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1107:11589:17295/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:1107:11589:17295/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1110:14096:95943/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:1110:14096:95943/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1110:15369:59046/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1110:15369:59046/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1111:8599:97362/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1111:8599:97362/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1113:10490:30873/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1113:10490:30873/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1113:12360:36316/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1113:12360:36316/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1113:4589:62685/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:1:1113:4589:62685/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1115:7288:99676/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1115:7288:99676/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1116:8136:52921/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1116:8136:52921/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1202:11809:66877/1	EQ|EQ|EQ	K01:2104=83/100M	K01:2104=83/100M	K01:2104=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1202:11809:66877/2	EQ|EQ|EQ	K01:2043=163/100M	K01:2043=163/100M	K01:2043=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1202:18844:98575/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1202:18844:98575/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1205:20782:28689/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1205:20782:28689/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1206:10108:83718/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:1206:10108:83718/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1212:17964:23344/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:1212:17964:23344/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1213:9111:56546/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1213:9111:56546/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1216:4380:98965/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1216:4380:98965/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1216:4493:85995/1	EQ|NE|NE	K01:2143=83/36S36M3I25M	K01:2143=83/36S36M3I25M	K01:2107=83/72M3I25M
HWI-1KL149:20:C1CU7ACXX:1:1216:4493:85995/2	EQ|EQ|EQ	K01:2043=163/100M	K01:2043=163/100M	K01:2043=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1216:8034:78319/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:1216:8034:78319/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:1316:14751:4679/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:1316:14751:4679/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2102:4725:60173/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:1:2102:4725:60173/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2104:19271:24502/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:1:2104:19271:24502/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2104:2016:81735/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2104:2016:81735/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2110:4445:72697/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:2110:4445:72697/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2111:2256:47748/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:2111:2256:47748/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2115:12497:79931/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:1:2115:12497:79931/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2115:17576:9737/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:2115:17576:9737/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2116:7977:30610/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2116:7977:30610/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2201:16984:100451/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2201:16984:100451/2	EQ|EQ|EQ	K01:2059=163/87M13S	K01:2059=163/87M13S	K01:2059=163/87M13S
HWI-1KL149:20:C1CU7ACXX:1:2203:19912:68616/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2203:19912:68616/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2204:13318:18341/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2204:13318:18341/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2206:10726:12303/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2206:10726:12303/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2206:11557:78671/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2206:11557:78671/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2211:12806:63973/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2211:12806:63973/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2212:17602:62052/1	EQ|EQ|EQ	K01:2104=83/100M	K01:2104=83/100M	K01:2104=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2212:17602:62052/2	EQ|EQ|EQ	K01:2043=163/100M	K01:2043=163/100M	K01:2043=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2212:19408:52552/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2212:19408:52552/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2303:8733:45438/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2303:8733:45438/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2304:9806:12935/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2304:9806:12935/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2305:12165:42334/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:1:2305:12165:42334/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2305:2388:67842/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2305:2388:67842/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2307:14199:91258/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2307:14199:91258/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2307:3121:93985/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2307:3121:93985/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2309:13907:13532/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2309:13907:13532/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2309:20396:57002/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:1:2309:20396:57002/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2312:11602:6630/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2312:11602:6630/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2313:11868:31327/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2313:11868:31327/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2313:9555:94108/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:1:2313:9555:94108/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:1:2315:19820:15046/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:1:2315:19820:15046/2	EQ|EQ|EQ	K01:2081=163/99M1S	K01:2081=163/99M1S	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1101:18362:28315/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1101:18362:28315/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1105:18846:45527/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1105:18846:45527/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1105:5659:65125/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1105:5659:65125/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1108:9609:39170/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1108:9609:39170/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1110:2262:8369/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1110:2262:8369/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1111:18496:5547/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1111:18496:5547/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1112:10132:23322/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1112:10132:23322/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1112:7260:56414/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1112:7260:56414/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1201:6906:82750/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:2:1201:6906:82750/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1202:16231:100362/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1202:16231:100362/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1213:12574:89489/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:2:1213:12574:89489/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1214:20898:3105/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1214:20898:3105/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:```20:C1CU7ACXX:2:1214:7035:46585/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1214:7035:46585/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1215:19107:31048/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1215:19107:31048/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1216:15500:73171/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:2:1216:15500:73171/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1216:6409:43952/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1216:6409:43952/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1301:16595:88662/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1301:16595:88662/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1306:12619:24138/1	EQ|NE|NE	K01:2143=83/36S36M3I25M	K01:2143=83/36S36M3I25M	K01:2107=83/72M3I25M
HWI-1KL149:20:C1CU7ACXX:2:1306:12619:24138/2	EQ|EQ|EQ	K01:2043=163/100M	K01:2043=163/100M	K01:2043=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1308:8618:21991/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:2:1308:8618:21991/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1309:15540:69632/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1309:15540:69632/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1314:9489:93274/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1314:9489:93274/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:1316:5692:99314/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:2:1316:5692:99314/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2108:14003:59876/1	EQ|EQ|EQ	K01:2136=83/44M3I53M	K01:2136=83/44M3I53M	K01:2136=83/44M3I53M
HWI-1KL149:20:C1CU7ACXX:2:2108:14003:59876/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2113:14713:81195/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2113:14713:81195/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2114:13002:89288/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2114:13002:89288/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2201:11170:94334/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:2:2201:11170:94334/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2202:4380:42920/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2202:4380:42920/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2213:14141:87844/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2213:14141:87844/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2213:15510:90052/1	EQ|EQ|EQ	K01:2122=83/57M3I23M1D17M	K01:2122=83/57M3I23M1D17M	K01:2122=83/57M3I23M1D17M
HWI-1KL149:20:C1CU7ACXX:2:2213:15510:90052/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2213:1759:92031/1	EQ|EQ|EQ	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M	K01:2123=83/56M3I41M
HWI-1KL149:20:C1CU7ACXX:2:2213:1759:92031/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2215:14279:12564/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2215:14279:12564/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2215:16028:26845/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2215:16028:26845/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2215:16774:44335/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2215:16774:44335/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2301:16926:82641/1	EQ|EQ|EQ	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M	K01:2136=83/43M3I54M
HWI-1KL149:20:C1CU7ACXX:2:2301:16926:82641/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2309:5711:82879/1	EQ|EQ|EQ	K01:2120=83/100M	K01:2120=83/100M	K01:2120=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2309:5711:82879/2	EQ|EQ|EQ	K01:1990=163/100M	K01:1990=163/100M	K01:1990=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2310:13450:61828/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2310:13450:61828/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2311:14082:99026/1	EQ|EQ|EQ	K01:2213=83/100M	K01:2213=83/100M	K01:2213=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2311:14082:99026/2	EQ|EQ|EQ	K01:2081=163/100M	K01:2081=163/100M	K01:2081=163/100M
HWI-1KL149:20:C1CU7ACXX:2:2315:4940:7934/1	EQ|EQ|EQ	K01:2133=83/100M	K01:2133=83/100M	K01:2133=83/100M
HWI-1KL149:20:C1CU7ACXX:2:2315:4940:7934/2	EQ|EQ|EQ	K01:2059=163/100M	K01:2059=163/100M	K01:2059=163/100M</h:pre>
```
END_DOC
*/

@Program(name="cmpbams",description="Compare two or more BAM files",
	keywords={"sam","bam","compare"}
	)
public class CompareBams  extends Launcher
	{
	private static final Logger LOG = Logger.build(CompareBams.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-Q","--mapq"},description="min MAPQ")
	private int min_mapq = 0 ;
	
	@Parameter(names={"-d","--distance"},description="distance tolerance between two alignments")
	private int distance_tolerance = 10 ;
	
	@Parameter(names={"-F","--sam"},description="use sam flag for comparaison")
	private boolean useSamFlag = false;
	
	@Parameter(names={"-c","--cigar"},description="use cigar String for comparaison")
	private boolean useCigar = false;
	
	@Parameter(names={"-r","--region"},description=IntervalParser.OPT_DESC)
	private String REGION = "";

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	

	
	
	private boolean samSequenceDictAreTheSame=true;
	private List<SAMSequenceDictionary> sequenceDictionaries=new ArrayList<SAMSequenceDictionary>();
	private PrintWriter out;
	
	private class MatchComparator
		implements Comparator<Match>
		{
		@Override
		public int compare(final Match m0, final Match m1)
			{
			int i=m0.readName.compareTo(m1.readName);
			if(i!=0) return i;
			i=m0.num_in_pair-m1.num_in_pair;
			if(i!=0) return i;
			//i= m0.bamIndex - m1.bamIndex;//NO ! (when comparing two Set<Match>)
			//if(i!=0) return i;
			i= compareTid(m0.bamIndex,m0.tid ,m1.bamIndex, m1.tid);
			if(i!=0) return i;
			i= m0.pos - m1.pos;
			if(i!=0) return i;
			i= m0.flag - m1.flag;
			if(i!=0) return i;
			i= m0.cigar.compareTo(m1.cigar);
			return 0;
			}
		}
	
	private class MatchOrderer
	implements Comparator<Match>
		{
		@Override
		public int compare(final Match m0, final Match m1)
			{
			int i=m0.readName.compareTo(m1.readName);
			if(i!=0) return i;
			i=m0.num_in_pair-m1.num_in_pair;
			return i;
			}
		}
	
	private class MatchCodec
		extends AbstractDataCodec<Match>
		{
		@Override
		public MatchCodec clone()
			{
			return new MatchCodec();
			}
		@Override
		public Match decode(final DataInputStream dis) throws IOException
			{
			final Match m=new Match();
			try {
				m.readName=dis.readUTF();
				}
			catch(IOException err)
				{
				return null;
				}
			m.bamIndex=dis.readInt();
			m.tid=dis.readInt();
			m.pos=dis.readInt();
			m.num_in_pair=dis.readInt();
			if(useSamFlag) m.flag=dis.readInt();
			if(useCigar) m.cigar=dis.readUTF();
			return m;
			}
		@Override
		public void encode(final DataOutputStream dos, final Match match)
				throws IOException
			{
			dos.writeUTF(match.readName);
			dos.writeInt(match.bamIndex);
			dos.writeInt(match.tid);
			dos.writeInt(match.pos);
			dos.writeInt(match.num_in_pair);
			if(useSamFlag) dos.writeInt(match.flag);
			if(useCigar) dos.writeUTF(match.cigar);
			}
		
		}
	
	private class Match
		{
		String readName;
		int num_in_pair=0;
		int tid=-1;
		int bamIndex=-1;
		int pos=-1;
		int flag=0;
		String cigar="";

		@Override
		public int hashCode()
			{
			int result = 1;
			result = 31 * result + num_in_pair;
			result = 31 * result + pos;
			result = 31 * result + tid;
			result = 31 * result + readName.hashCode();
			result = 31 * result + bamIndex;
			result = 31 * result + cigar.hashCode();
			return result;
			}
		@Override
		public boolean equals(final Object obj)
			{
			if (this == obj) { return true; }
			if (obj == null) { return false; }
			final Match other = (Match) obj;
			if (num_in_pair != other.num_in_pair) { return false; }
			if (compareTid(this.bamIndex,tid,other.bamIndex,other.tid)!=0) { return false; }
			if(tid==-1) return true;
			if (pos != other.pos) { return false; }
			if (bamIndex != other.bamIndex) { return false; }
			if(!readName.equals(other.readName)) return false;
			return true;
			}
		}
	
	private String norm(String s1)
		{
		if(s1.startsWith("chr")) s1=s1.substring(3);
		if(s1.startsWith("0")) s1=s1.substring(1);
		if(s1.equals("MT")) s1="M";
		return s1;
		}
	
	private int compare(final String s1,final String s2)
		{
		return norm(s1).compareToIgnoreCase(norm(s2));
		}
	
	private int compareTid(int file_id1,int tid1,int file_id2,int tid2)
		{
		if(samSequenceDictAreTheSame) return tid1-tid2;
		if(tid1==-1)
			{
			return tid2==-1?0:-1;
			}
		if(tid2==-1)
			{
			return 1;
			}
		final String chrom1=this.sequenceDictionaries.get(file_id1).getSequence(tid1).getSequenceName();
		final String chrom2=this.sequenceDictionaries.get(file_id2).getSequence(tid2).getSequenceName();
		if(chrom1==null)
			{
			return chrom2==null?0:-1;
			}
		if(chrom2==null) return 1;
		return compare(chrom1,chrom2);
		}
	
	private void print(final Set<Match> set,final SAMSequenceDictionary dict)
		{
		boolean first=true;
		for(final Match m:set)
			{
			if(!first)this.out.print(',');
			first=false;
			if(m.tid<0){ this.out.print("unmapped"); continue;}
			final SAMSequenceRecord ssr=(dict==null?null:dict.getSequence(m.tid));
			String seqName=(ssr==null?null:ssr.getSequenceName());
			if(seqName==null) seqName="tid"+m.tid;
			this.out.print(String.valueOf(seqName+":"+(m.pos)));
			if(this.useSamFlag) this.out.print("="+m.flag);
			if(this.useCigar) this.out.print("/"+m.cigar);
			}
		if(first) this.out.print("(empty)");
		}
	
	
	
	
	
	private final List<File> IN=new ArrayList<File>();
	
    private boolean same(final Set<Match> set1,final Set<Match> set2)
    	{
    	for(final Match m0:set1)
    		{
    		for(final Match m1:set2)
	    		{
    			int i=m0.readName.compareTo(m1.readName);
    			if(i!=0) continue;
    			i=m0.num_in_pair-m1.num_in_pair;
    			if(i!=0) continue;
    			//i= m0.bamIndex - m1.bamIndex;//NO ! (when comparing two Set<Match>)
    			//if(i!=0) return i;
    			i= m0.tid - m1.tid;
    			if(i!=0) continue;
    			i= Math.abs(m0.pos - m1.pos);
    			if(i>this.distance_tolerance) continue;
    			i= m0.flag - m1.flag;
    			if(i!=0) continue;
    			i= m0.cigar.compareTo(m1.cigar);
    			if(i!=0) continue;
    			return true;
	    		}
    		}
    	return false;
    	}
    
    @Override
    public int doWork(final List<String> args) {
    	this.IN.addAll(args.stream().map(S->new File(S)).collect(Collectors.toList()));
   		SortingCollection<Match> database = null;
		SamReader samFileReader=null;
		CloseableIterator<Match> iter=null;
		try
			{
			if(this.IN.size() <2)
				{
				LOG.error("Need more bams please");
				return -1;
				}
			
			database = SortingCollection.newInstance(
					Match.class,
					new MatchCodec(),
					new MatchOrderer(),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			this.samSequenceDictAreTheSame=true;
			database.setDestructiveIteration(true);
	
			
			for(int currentSamFileIndex=0;
					currentSamFileIndex<this.IN.size();
					currentSamFileIndex++ )
				{
				final File samFile=this.IN.get(currentSamFileIndex);
				LOG.info("Opening "+samFile);
				samFileReader= super.createSamReaderFactory().open(samFile);
				final SAMSequenceDictionary dict=samFileReader.getFileHeader().getSequenceDictionary();
				if(dict==null || dict.isEmpty())
					{
					LOG.error("Empty Dict  in "+samFile);
					return -1;
					}
				
				if(!this.sequenceDictionaries.isEmpty() &&
					!SequenceUtil.areSequenceDictionariesEqual(this.sequenceDictionaries.get(0), dict))
					{
					this.samSequenceDictAreTheSame=false;
					LOG.warn("FOOL !! THE SEQUENCE DICTIONARIES ARE **NOT** THE SAME. I will try to compare anyway but it will be slower.");
					}
				this.sequenceDictionaries.add(dict);
				
				
				final Optional<Interval> interval;
				if(REGION!=null && !REGION.trim().isEmpty())
					{
					final IntervalParser dix = new IntervalParser(dict);
					interval = Optional.ofNullable(dix.parse(REGION));
					
					if(!interval.isPresent())
						{
						LOG.error("Cannot parse "+REGION+" (bad syntax or not in dictionary)");
						return -1;
						}
					}
				else
					{
					interval = Optional.empty();
					}
				
				
				SAMRecordIterator it=null;
				if(!interval.isPresent())
					{
					it=samFileReader.iterator();
					}
				else
					{
					it=samFileReader.queryOverlapping(
							interval.get().getContig(),
							interval.get().getStart(),
							interval.get().getEnd()
							);
					}
				final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
				while(it.hasNext() )
					{
					final SAMRecord rec=progress.watch(it.next());
					if(!rec.getReadUnmappedFlag())
						{
						if(rec.getMappingQuality() < this.min_mapq) continue;
						if(rec.isSecondaryOrSupplementary()) continue;
						}
					final Match m=new Match();
					if(rec.getReadPairedFlag())
						{
						m.num_in_pair=(rec.getFirstOfPairFlag()?1:2);
						}
					else
						{
						m.num_in_pair=0;
						}
					m.readName=rec.getReadName();
					m.bamIndex=currentSamFileIndex;
					m.flag=rec.getFlags();
					m.cigar=rec.getCigarString();
					if(m.cigar==null ) m.cigar="";
					if(rec.getReadUnmappedFlag())
						{
						m.tid=-1;
						m.pos=-1;
						}
					else
						{
						m.tid=rec.getReferenceIndex();
						m.pos=rec.getAlignmentStart();
						}
					database.add(m);
					}
				it.close();
				samFileReader.close();
				samFileReader=null;
				LOG.info("Close "+samFile);
				}
			database.doneAdding();
			LOG.info("Writing results....");
			
			this.out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			//compute the differences for each read
			this.out.print("#READ-Name\t");
			for(int x=0;x<this.IN.size();++x)
				{
				for(int y=x+1;y<this.IN.size();++y)
					{
					if(!(x==0 && y==1)) this.out.print("|");
					this.out.print(IN.get(x));
					this.out.print(" ");
					this.out.print(IN.get(y));
					}
				}
			for(int x=0;x<this.IN.size();++x)
				{
				this.out.print("\t"+IN.get(x));
				}
			this.out.println();
			
			/* create an array of set<Match> */
			final MatchComparator match_comparator=new MatchComparator();
			final List<Set<Match>> matches=new ArrayList<Set<CompareBams.Match>>(this.IN.size());
			while(matches.size() < this.IN.size())
				{
				matches.add(new TreeSet<CompareBams.Match>(match_comparator));
				}
			
			iter = database.iterator();
			String currReadName=null;
			int curr_num_in_pair=-1;
			for(;;)
				{
				Match nextMatch = null;
				if(iter.hasNext())
					{
					nextMatch = iter.next();
					}
				if(nextMatch==null ||
					(currReadName!=null && !currReadName.equals(nextMatch.readName)) ||
					(curr_num_in_pair!=-1 && curr_num_in_pair!=nextMatch.num_in_pair))
					{
					if(currReadName!=null)
						{
						this.out.print(currReadName);
						if(curr_num_in_pair>0)
							{
							this.out.print("/");
							this.out.print(curr_num_in_pair);
							}
						this.out.print("\t");
						
						
						for(int x=0;x<this.IN.size();++x)
							{
							final Set<Match> first=matches.get(x);
							for(int y=x+1;y<this.IN.size();++y)
								{
								if(!(x==0 && y==1)) this.out.print("|");
								Set<Match> second=matches.get(y);
								if(same(first,second))
									{
									this.out.print("EQ");
									}
								else
									{
									this.out.print("NE");
									}
								}
							}
	
						for(int x=0;x<this.IN.size();++x)
							{
							this.out.print("\t");
							print(matches.get(x),sequenceDictionaries.get(x));
							}
						
						this.out.println();
						}
					if(nextMatch==null) break;
					for(Set<Match> set:matches) set.clear();
					}
				currReadName=nextMatch.readName;
				curr_num_in_pair=nextMatch.num_in_pair;
				matches.get(nextMatch.bamIndex).add(nextMatch);
				if(this.out.checkError()) break;
				}
			
			iter.close();
			this.out.flush();
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(database!=null) database.cleanup();
			CloserUtil.close(samFileReader);
			CloserUtil.close(this.out);this.out=null;
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		new CompareBams().instanceMainWithExit(args);
		}
}
