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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamJsonWriterFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.util.RuntimeIOException;
/**
BEGIN_DOC

## Example

```
$ java -jar dist/sam2json.jar src/test/resources/toy.bam | python -m json.tool
[
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            },
            {
                "name": "XX",
                "value": [
                    12561,
                    2,
                    20,
                    112
                ]
            }
        ],
        "cigar": "8M4I4M1D3M",
        "flag": 163,
        "len": 39,
        "mapq": 30,
        "matepos": 37,
        "materef": "ref",
        "name": "r001",
        "pos": 7,
        "qualities": "*",
        "ref": "ref",
        "sequence": "TTAGATAAAGAGGATACTG"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "1S2I6M1P1I1P1I4M2I",
        "flag": 0,
        "mapq": 30,
        "name": "r002",
        "pos": 9,
        "qualities": "*",
        "ref": "ref",
        "sequence": "AAAAGATAAGGGATAAA"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "5H6M",
        "flag": 0,
        "mapq": 30,
        "name": "r003",
        "pos": 9,
        "qualities": "*",
        "ref": "ref",
        "sequence": "AGCTAA"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "6M14N1I5M",
        "flag": 0,
        "mapq": 30,
        "name": "r004",
        "pos": 16,
        "qualities": "*",
        "ref": "ref",
        "sequence": "ATAGCTCTCAGC"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "6H5M",
        "flag": 16,
        "mapq": 30,
        "name": "r003",
        "pos": 29,
        "qualities": "*",
        "ref": "ref",
        "sequence": "TAGGC"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "9M",
        "flag": 83,
        "len": -39,
        "mapq": 30,
        "matepos": 7,
        "materef": "ref",
        "name": "r001",
        "pos": 37,
        "qualities": "*",
        "ref": "ref",
        "sequence": "CAGCGCCAT"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "20M",
        "flag": 0,
        "mapq": 30,
        "name": "x1",
        "pos": 1,
        "qualities": "????????????????????",
        "ref": "ref2",
        "sequence": "AGGTTTTATAAAACAAATAA"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "21M",
        "flag": 0,
        "mapq": 30,
        "name": "x2",
        "pos": 2,
        "qualities": "?????????????????????",
        "ref": "ref2",
        "sequence": "GGTTTTATAAAACAAATAATT"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "9M4I13M",
        "flag": 0,
        "mapq": 30,
        "name": "x3",
        "pos": 6,
        "qualities": "??????????????????????????",
        "ref": "ref2",
        "sequence": "TTATAAAACAAATAATTAAGTCTACA"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "25M",
        "flag": 0,
        "mapq": 30,
        "name": "x4",
        "pos": 10,
        "qualities": "?????????????????????????",
        "ref": "ref2",
        "sequence": "CAAATAATTAAGTCTACAGAGCAAC"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "24M",
        "flag": 0,
        "mapq": 30,
        "name": "x5",
        "pos": 12,
        "qualities": "????????????????????????",
        "ref": "ref2",
        "sequence": "AATAATTAAGTCTACAGAGCAACT"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "23M",
        "flag": 0,
        "mapq": 30,
        "name": "x6",
        "pos": 14,
        "qualities": "???????????????????????",
        "ref": "ref2",
        "sequence": "TAATTAAGTCTACAGAGCAACTA"
    }
]
```

END_DOC
 */
@Program(name="sam2json",
	keywords={"sam","bam","json"},
	description="Convert a SAM input to JSON",
	modificationDate="20210315",
	creationDate="20210402")
public class SamToJson extends OnePassBamLauncher {
	private static final Logger LOG = Logger.build(SamToJson.class).make();
	
	@Parameter(names={"-H","--header"},description="don't print SAM HEADER")
	private boolean print_header = false;

	@Parameter(names={"-name","--name"},description="do not print read name")
	private boolean disable_readName = false;

	@Parameter(names={"-atts","--atts"},description="do not print attributes")
	private boolean disable_atts = false;

	@Parameter(names={"-flag","--flag"},description="expand SAm Flags")
	private boolean expflag = false;

	@Parameter(names={"-cigar","--cigar"},description="expand cigar")
	private boolean excigar = false;

	@Override
	protected Logger getLogger() {
		return LOG;
		}

	@Override
	protected SAMFileWriter openSamFileWriter(SAMFileHeader headerIn)
		{
		final SamJsonWriterFactory factory =SamJsonWriterFactory.newInstance().
				printHeader(this.print_header).
				printReadName(!this.disable_readName).
				printAttributes(!this.disable_atts).
				expandFlag(this.expflag).
				expandCigar(this.excigar)
				;
		PrintWriter out;
		try {
			out= super.openPathOrStdoutAsPrintWriter(super.outputFile);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		return factory.open(super.createOutputHeader(headerIn), out);
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamToJson().instanceMainWithExit(args);

	}

}
