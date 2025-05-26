
/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
* 2015 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.nio.file.Path;
import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

/**
BEGIN_DOC

## See also

* https://twitter.com/sjackman/status/584418230791340032

## Example

```
$ make bamtile && \
  curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110915_CEUtrio_b37_decoy_alignment/CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam" |\
  java -jar dist/jvarkit.jar bamtile
(...)
B00EGABXX110201:5:25:14956:33658	121	1	9994	23	101M	=	9994	0	CTTCCGATCTCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC	AC:C;BA<CACCDCD=CACBD9?AE?D?CDECDDDDECDDDDECDDDDEDD?DDECDCDDDCECDDECECDDEDECDDDBDADDD?DBDDD@DACCDAC@A	X0:i:1	X1:i:1	XA:Z:X,+155260188,101M,5;	MD:Z:0T0G6A0A91	RG:Z:B00EG.5	XG:i:0	AM:i:0	NM:i:4	SM:i:23	XM:i:4	XN:i:7	XO:i:0BQ:Z:\^U^V]\W^\PPQPQJPNPOQFLNRLQLPQRPQQQQRPQQQQRPQQQQRQQLQQRPQPQQQPRPQQRPRPQQRQRPQQQOQNQQQLQOQQQMQNPPQNPMN	OQ:Z:EE:EEEC;EACEEEE7DBCCD/<>@>D6ADDEFBFEDCE>EFECE@FFFFF9EEDDG@GFGFHEHHHFHFHHHHHHGGGFGCGGG@HFHHHEHFHHHGHHH	XT:A:U
B06PYABXX110322:7:2203:15946:92062	99	1	10091	1	87M1I3M1I7M2S	=	10378	371	TAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCAAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCCTAAA	AACACCDCEBCCDCEBDBCECDDEBDCDDEDEB@DDCECCDEDFBDCBEDB<CECE@CCCECCBAAEBDBDCD@DBCCE7D@@@C?DBACC76<A4>D###	X0:i:1	X1:i:3XA:Z:1,+10091,87M1I3M1I7M,5;21,-48119790,5M1I33M1I59M,5;1,-249240224,5M1I27M1D6M1I59M,5;	XC:i:99	MD:Z:47T49	RG:Z:B06PY.7	XG:i:2	AM:i:1	NM:i:3	SM:i:1	XM:i:2	XO:i:1	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\^[@QU@@@@@@@	MQ:i:1	OQ:Z:HHHHGHGHHHHHHGHHGEHHGHHHEFHHHHHFF@HFFHHFHDHHDHFBFEF<FEEFBDFEEFBFBBDCEEDC=?EDDEE5DAAAA=EEBF?55:C5=@###	XT:A:U
B00EGABXX110201:6:47:12319:33409	99	1	10185	29	63M1D38M	=	10363	267	CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTACCCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTA	ABADBCCDBEBCCDCECDDECECDDECFCDDECECDDE@ECCDDECBDDDECECCDF=E@CDFCA<CBFAE?CC@ECEBCDE?EBC=D>D@BAA5:579;6	X0:i:3	X1:i:4MD:Z:43T0A18^A8A29	RG:Z:B00EG.6	XG:i:1	AM:i:0	NM:i:4	SM:i:0	XM:i:3	XO:i:1	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Z]__`]@@@@@@@@@@@@@@@@@^[@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:37	OQ:Z:HHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHHHHAHHFHHFHFFHHCHFHEEB=GDDECB=@CCAACACCBEFEFFHH?HHF@A?EDDC;57:;<46	XT:A:R
B00EGABXX110201:5:3:9889:46988	163	1	10285	29	49M1D52M	=	10398	177	TAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTACCCC	>AA@ABCBCACCCBCBBCCBD9BCB?EABCCCDBBB=ACACBCBDBCC6C>@BB=BC?B9/A3<98:898679:6B<=B7>AA7@39858795281?588,	X0:i:1	X1:i:0	MD:Z:6C35T6^T5T42A3	RG:Z:B00EG.5	XG:i:1	AM:i:29	NM:i:5	SM:i:29	XM:i:4	XO:i:1	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@GKMA@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:29	OQ:Z:HHEHHHHHFHHHHFFHFHHFFBFGFCHFFHHHFHFFBDEFFFFEFFFF<F>EEEBEEEE>*D5C@>6<B?=2=<>EA:E9DDD2D5A>=5;;@:@+C@@?3	XT:A:U
C00D5ABXX110314:1:1104:14854:27916	147	1	10368	29	1S18M2D2M1D6M1D52M1I20M1S	=	10195	-274	ACTACCCCTAACCCTAACCAACCCTAACCCTAACCCTAACACTAACCCTAACCCTCACCCTCACCCTAACCCTTAACCCTTAACCCTAACCCTCACCCTCA	#@;@,=68.6070316>;@>:BAD@D4BBDCB=ECDEFDECEEFDDDAEEDDDDDD:CCDDDCCADDEDBDDCDDC@DDEDDBDDCCD@CCBCB7ABBCA?	MD:Z:3A14^CT2^C6^C13C14A5A10C19A6	RG:Z:C00D5.1	XG:i:5	AM:i:29	NM:i:11	SM:i:29	XM:i:6	XO:i:4	BQ:Z:@[V[GXQSIQKRKNLQYV[YU]\_[_O]]_^]X`^_`a_`^``a___\``______U^^___^^\__`_]__^__^[__@__]__^^_[^^]^]R\]]^\@	MQ:i:29OQ:Z:#@:?0?65-3'83113<?>:7DDD@C0DDDBB<GEEGGFGDGGGFGHEHHFFFFCF:EEEDHEFCIFGEDGGBGEFCGGFGFCGGFEGAGGEDE7EGHHHH	XT:A:M
B00EGABXX110201:6:27:4809:82569	163	1	10467	60	101M	=	10517	150	CTCGCGGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACAACGCAG	>BA9B:@BAA?CD@BCB6:C'=<5=AB92B66A=;>DCA@BDA8<BBC+C@EDCDDE:C:?A@BA@@BA7B3?B@=8=-=@=@;>A7;ADCB4B:C?:C;7	X0:i:1	X1:i:0	MD:Z:101	RG:Z:B00EG.6	XG:i:0	AM:i:37	NM:i:0	SM:i:37	XM:i:0	XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:60	OQ:Z:HEGHHGCGGGCHGDGFF=DF.AB=BEEA6E=<C@:BEEBEEEB9;DDD*DEFFFFGG?FECD>DDBADC4E=CEEB@?.@CB@;@C=ACGEE6G?EEDG@5	XT:A:U
C00D5ABXX110314:6:1105:12956:64901	147	1	10525	15	101M	=	10316	-309	CGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCGGCGC	9796<37;:4:?2=@6;8>BB8ADDC>@CDCECA;DA6;BBCAEAD<CDDD?<D;A@CBCDDD?8CC>EDCA<BCDBDB?DBCB;@BDACAC;BB:BA9A>	X0:i:1	X1:i:1	XA:Z:12,+94998,101M,3;	MD:Z:97T0T2	RG:Z:C00D5.6	XG:i:0	AM:i:15	NM:i:2	SM:i:15	XM:i:2	XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@D@DA	MQ:i:15	OQ:Z:A6;181157.:=+>=1;?<A>2<EEG;@CCFDDA7DB3ECBBBH?HHFHFF8EECBABBCBEDB?EHBEHFF:CFHEFFBDDECEBFFBFAHHHHHHHHHH	XT:A:U
B00EGABXX110201:6:3:7973:149498	83	1	10740	0	66S35M	=	10473	-301	GGCCCCCCCCTCGCGGTACTTGCTGGGTCTGTTGTGAGGAGAACGCAACTCCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGAG	###################################################################=E<E?1<A8BC;>9A?8BB;6CCBB:BBCAA?@>	X0:i:39	XC:i:35	MD:Z:35	RG:Z:B00EG.6	XG:i:0	AM:i:0	NM:i:0	SM:i:0	XM:i:0	XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@>X`W`ZLW\S]^VYT\ZS]]VQ^^]]U]]^\\Z[Y	MQ:i:37	OQ:Z:###################################################################?E<DA2EC?EED@BDB?EE@7FFFFEFFFFDFFF	XT:A:R
B00EGABXX110201:4:24:10501:106133	99	1	10760	0	35M66S	=	11009	331	CGCAGGCGCAGAGAGGCGCGCCGCGGCGGCGCAAGCGCAGAGACACATGCTCCCGGGACCAAGGGTCGGGGGCGGGGGGCGCGGAGCGGCCGCCCACCACC	8437>8794;<9>*B<7;###################################################################################	X0:i:50	XC:i:35	MD:Z:25C7G1	RG:Z:B00EG.4	XG:i:0	AM:i:0	NM:i:2	SM:i:0	XM:i:2	XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:0	OQ:Z:<66:96:<6866:)@<8@###################################################################################	XT:A:R
B06PYABXX110322:6:1208:11671:56074	163	1	10782	0	35M66S	=	10883	186	GCGCCGGCGCAGGCGCAGACACATGCTAGCGCGTCCGGGGGGGAGGCGTGGCGCAGGAGCAGAGAGGCGCGGAGCGCCCGCGGAGGCGCAGAGACAGATGC	<A9A:7CB;CBCCB:CBCDBCBB0C?23?########################################################################	X0:i:2	X1:i:1	XA:Z:1,+10782,101M,0;hs37d5,+9321764,101M,1;	XC:i:35	MD:Z:35	RG:Z:B06PY.6	XG:i:0	AM:i:0	NM:i:0	SM:i:0	XM:i:0XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:0	OQ:Z:BEEE=BGFGGHEHFFGDDGGGGF-D?.5?########################################################################	XT:A:R
B00EGABXX110201:5:44:9232:74422	163	1	10824	0	101M	=	11008	258	AGGCGTGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAGAGGC	>BAA9@BCB:CBDCA:BCDCDBCC=:C8@B81:7B8?4/<<A>+29:?<@=8::=BDCBB9<6<;;BA;=C?=9?C7A94/182&($'3587-1?>@6*=#	X0:i:2	X1:i:1	XA:Z:15,-102520246,101M,0;hs37d5,+9321807,101M,1;	MD:Z:101	RG:Z:B00EG.5	XG:i:0	AM:i:0	NM:i:0	SM:i:0	XM:i:0XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:0	OQ:Z:HFFGFEFHGHHGHHEFFHHFHEEFBCFBDEB7C<EAB92AABA17>>?>><?>BAEFFCEBB?AE9EDA=EBA7@E:C<93*76,.&'19>@33>AA:(A#	XT:A:R
C00D5ABXX110314:5:2102:8433:99141	99	1	10923	0	87M14S	=	10968	129	GCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAGACGCAAGCCTACGGGCGGGGGGTGGGG	BA:CCBC;C<DD<DD<DDEDC<DEEEDEBDCCDCDEDEC<D<CDDDDDDD?DD>C@B;>D@C<?4C>+:@-A/B..65,9=5850A###############	X0:i:2	X1:i:1	XA:Z:1,+10923,101M,0;hs37d5,+9321906,101M,1;	XC:i:87	MD:Z:87	RG:Z:C00D5.5	XG:i:0	AM:i:0	NM:i:0	SM:i:0	XM:i:0XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:0	OQ:Z:HHHHHHHHHGGHHGHHHFHHFHHGHFFEFFGCHFFEFHHHHHBHHEEEFF:DE9>A>A7AAAA;.??,<;';*:*.63(3637/.A###############	XT:A:R
B00EGABXX110201:4:24:10501:106133	147	1	11009	0	18S83M	=	10760	-331	CAGAGGCCGCAGCGCAACGGGCGGGGGTTGGGGGGGCGTGTGTGGCAGGAGCAAAGTCGCACGGCGCCGGGCTGGGGCGGGGGGAGGGTGGCGCCGTGCAC	###################95&<</<$17599%A?>76?B?>A6A=DBBBC@=>>9>4>>>8B:;CD6BBC><<?<A8<B:BBA@BCBBCC:BB:@A@A@>	X0:i:2	X1:i:0	XA:Z:15,+102520079,101M,1;	XC:i:83	MD:Z:25T57	RG:Z:B00EG.4	XG:i:0	AM:i:0	NM:i:1	SM:i:0	XM:i:1	XO:i:0BQ:Z:@@@@@@@@@@@@@@@@@@>TPAWWJW?LRPTT@\ZYRQZ]ZY\Q\X_]]]^[XYYTYOYYYS]UV^_Q]]^YWWZW\SW]U]]\[]^]]^^U]]U[\[\[Y	MQ:i:0	OQ:Z:###################@;)AA6A./;=?>(DCB@@BFCDD=D=EEEBE>@?:A@;BB=@E?EFG?FFF=?BCBDBBE>FFD?FHHGHGFGIEGGEHHG	XT:A:R
B06PYABXX110322:7:2107:20624:6619	99	1	11082	0	101M	=	11276	288	CGCCGTGCACGCGCAGAAACTCACGTCACGGTGGCGCGGCGCAGAGACGGGGAGAACCTCAGTAATCCGAAAAGCCGGGATCGACCGCCCCTTGCTTGCAG	@9AB;BCCCB<C<CCDDEEBCCDC;ADCC=CBDCC=C<DC;DCDADB?<DC#B@=B?AC>>@8@>?CB;2CBCCC>8??A<A:AB8::C@26>@=@?B<B#	X0:i:3	X1:i:0	XA:Z:1,+11082,101M,1;hs37d5,+9322064,101M,1;	MD:Z:51T49	RG:Z:B06PY.7	XG:i:0	AM:i:0	NM:i:1	SM:i:0	XM:i:1XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:0	OQ:Z:GGGGGGGGGGGGGGGGGGGGCFGGEAGFGGECGDFGCGFCCDDF>EABFED(A;7A?A=;9;.>:9DAC.A@A@>??AA?4AA>D7B7D@1.:::9<>9B#	XT:A:R
B00EGABXX110201:5:3:10660:82271	83	1	11181	0	101M	=	11111	-170	AGCCGGGCACTACAGGACCCGCTTGCTCACGGTGCTGTGCCAGGGCGCCCCCTGCTGGCGACTAGGGCAACTGCAGGGCTCTCTTGCTTAGAGTGGTGGCC	#?A7A?;@@@?5??:4AA'<ADCCDDDDE<DBCDECBCCDDECDD<9DDDDEDDECDD=DCECEDDDDECECDDEDDDEDEDEECDEECDCDBBCBBBBBA	X0:i:7	X1:i:0	MD:Z:101	RG:Z:B00EG.5	XG:i:0	AM:i:0	NM:i:0	SM:i:0	XM:i:0	XO:i:0	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	MQ:i:0	OQ:Z:#BC;CB??<<?.:7=5<B,@BFBFFFHHBHHHHHHHHBDFFFHHFF9FFGFHFHHHHHGHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHH	XT:A:R
```

print name and prev.end-current.start with awk:

```
curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110915_CEUtrio_b37_decoy_alignment/CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam" |\
java -jar dist/jvarkit.jar bamtile | grep -v '^@' |\
awk -F ' ' 'BEGIN{prev=-1;} {printf("%s %d\n",$1,int($4)-prev);prev=int($4)+length($10);}' 
B00EGABXX110201:5:25:14956:33658 9995
B06PYABXX110322:7:2203:15946:92062 -4
B00EGABXX110201:6:47:12319:33409 -7
B00EGABXX110201:5:3:9889:46988 -1
C00D5ABXX110314:1:1104:14854:27916 -18
B00EGABXX110201:6:27:4809:82569 -2
C00D5ABXX110314:6:1105:12956:64901 -43
B00EGABXX110201:6:3:7973:149498 114
B00EGABXX110201:4:24:10501:106133 -81
B06PYABXX110322:6:1208:11671:56074 -79
B00EGABXX110201:5:44:9232:74422 -59
C00D5ABXX110314:5:2102:8433:99141 -2
B00EGABXX110201:4:24:10501:106133 -15
B06PYABXX110322:7:2107:20624:6619 -28
B00EGABXX110201:5:3:10660:82271 -2
B06PYABXX110322:8:2208:13299:179018 -4
B00EGABXX110201:6:1:18827:122731 -7
B00EGABXX110201:5:68:12590:13115 -1
C00D5ABXX110314:5:2201:6237:9231 -1
B00EGABXX110201:5:8:4957:37361 -1
C00D5ABXX110314:5:1208:12083:200469 -1
B00EGABXX110201:5:7:6075:93915 -1
B00EGABXX110201:5:27:16780:84550 -4
B06PYABXX110322:6:1107:4742:50543 -11
B06PYABXX110322:7:1205:12471:106838 -1
B00EGABXX110201:5:26:7748:160165 -2
B00EGABXX110201:4:27:18455:19154 -17
B06PYABXX110322:5:2101:14995:168978 -11
B00EGABXX110201:5:61:10922:47487 -4
C00D5ABXX110314:2:1201:10252:100348 -4
B06PYABXX110322:8:1108:2771:165330 -7
B06PYABXX110322:8:2107:4542:82481 -1
B06PYABXX110322:5:2107:16237:128403 -8
B00EGABXX110201:6:2:7795:198129 -12
B06PYABXX110322:8:1106:17730:27458 -5
B06PYABXX110322:8:1204:15694:112071 -1
C00D5ABXX110314:6:1208:8777:72655 -1
B06PYABXX110322:5:2201:2780:172842 -2
B06PYABXX110322:5:2203:9768:100733 -2
B06PYABXX110322:5:2106:4346:156704 -4
B06PYABXX110322:5:2206:18887:48252 -3
B06PYABXX110322:6:2107:10460:183803 -2
B06PYABXX110322:5:2204:6403:169005 -11
B00EGABXX110201:5:5:7775:24618 -14
C00D5ABXX110314:2:2104:4859:131692 -1
B06PYABXX110322:6:2108:18994:79163 -2
B06PYABXX110322:6:1202:11434:192811 -3
B06PYABXX110322:5:2207:5218:123898 -1
B06PYABXX110322:6:2108:10847:149590 -9
B06PYABXX110322:8:1204:19832:68914 -3
```


END_DOC
*/
@Program(
		name="bamtile",
		description="Answer to @sjackman : Is there a bedtools command to determine a minimal tiling path? A minimal set of features that cover a maximum of the target.",
		keywords={"bam","sam","tile"},
		biostars=287915,
		modificationDate="20191010",
		jvarkit_amalgamion = true
		)
public class BamTile
	extends Launcher
	{
	private static final Logger LOG = Logger.of(BamTile.class);

	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description="For Reading Cram. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx = null;
	@Parameter(names={"-e","--exclude"},description="[20171206]"+SamRecordJEXLFilter.FILTER_DESCRIPTION)
	private SamRecordFilter filterOut = SamRecordJEXLFilter.buildAcceptAll();
	@Parameter(names={"-n","--no-overlap"},description="[20171206]No overlap, just the read close to each other")
	private boolean no_overlap  = false;

	public BamTile()
		{
		}

	@Override
	public int doWork(final List<String> args)
		{					
		SAMRecordIterator iter=null;
		
		SamReader sfr=null;
		SAMFileWriter sfw =null;
		try
			{			
			sfr = openSamReader(oneFileOrNull(args));
			final SAMFileHeader header1=sfr.getFileHeader();
			if(header1==null)
				{
				LOG.error("File header missing");
				return -1;
				}
			
			if(header1.getSortOrder()!=SAMFileHeader.SortOrder.coordinate)
				{
				LOG.error("File header not sorted on coordinate");
				return -1;
				}
			
			final SAMFileHeader header2=header1.clone();
			JVarkitVersion.getInstance().addMetaData(this, header2);			
			sfw =  this.writingBamArgs.setReferencePath(this.faidx).openSamWriter(this.outputFile,header2, true);
			
			
			final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.newInstance().dictionary(header1).logger(LOG).build();
			iter=sfr.iterator();
			final LinkedList<SAMRecord> buffer=new LinkedList<>();
			for(;;)
				{
				SAMRecord rec=null;
				if( iter.hasNext())
					{
					rec= progress.apply(iter.next());
					if(rec.getReadUnmappedFlag()) continue;
					if(this.filterOut.filterOut(rec)) continue;
					if(!buffer.isEmpty())
						{
						final SAMRecord last= buffer.getLast();
						if(this.no_overlap)
							{
							if(last.getReferenceIndex()==rec.getReferenceIndex()  && 
							   last.getEnd()>=rec.getStart())
								{
								continue;
								}
							}
						else
							{
							if( last.getReferenceIndex()==rec.getReferenceIndex() &&
								last.getAlignmentStart()<=rec.getAlignmentStart() &&
								last.getAlignmentEnd()>=rec.getAlignmentEnd())
								{
								continue;
								}
							}
						}
					}
				if(rec==null || (!buffer.isEmpty() && buffer.getLast().getReferenceIndex()!=rec.getReferenceIndex()))
					{
					while(!buffer.isEmpty())
						{
						sfw.addAlignment(buffer.removeFirst());
						}
					if(rec==null) break;
					}
				buffer.add(rec);
				
				
				if(!this.no_overlap && buffer.size()>2)
					{
					final int index = buffer.size();
					final SAMRecord prev =  buffer.get(index-3);
					final SAMRecord curr =  buffer.get(index-2);
					final SAMRecord next =  buffer.get(index-1);
					
					if( prev.getAlignmentEnd() >= next.getAlignmentStart() ||
						curr.getAlignmentEnd() <= prev.getAlignmentEnd())
						{
						buffer.remove(index-2);
						}
					else if(curr.getAlignmentStart() == prev.getAlignmentStart() &&
							curr.getAlignmentEnd() > prev.getAlignmentEnd()
							)
						{
						buffer.remove(index-3);
						}
					}
				while(buffer.size()>3)
					{
					sfw.addAlignment(buffer.removeFirst());
					}
				
				}
			progress.close();
			sfw.close();sfw=null;
			iter.close();iter=null;
			sfr.close();sfr=null;
			return  RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
				
	public static void main(final String[] args)
		{
		new BamTile().instanceMainWithExit(args);
		}
	}
