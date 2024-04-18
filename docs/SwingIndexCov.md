# SwingIndexCov

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

indexcov visualization


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar swingindexcov  [options] Files

Usage: swingindexcov [options] Files
  Options:
    --gtf, --gff, --gff3
      GFF3 file indexed with tabix to plot the genes.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --helper
      For expert users only. java archive implenting Helper. Syntax 
      "path/to/helper.jar package.helper.implementation.Name"
    -R, --reference
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --version
      print version and exit

```


## Keywords

 * cnv
 * duplication
 * deletion
 * sv



## Creation Date

2020511

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/SwingIndexCov.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/SwingIndexCov.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **swingindexcov** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Examples:

```
java -jar dist/swingindexcov.jar indexcov.bed.gz
```
or 
```
java -jar dist/swingindexcov.jar -R reference.fa indexcov.bed.gz --gff indexed.gff3.gz --helper "/TMP/myhelper.jar MyHelper"
```

## Screenshots


![https://pbs.twimg.com/media/FSkLDJdXwAEmZTO?format=png&name=small](https://pbs.twimg.com/media/FSkLDJdXwAEmZTO?format=png&name=small
)

![https://pbs.twimg.com/media/FSkLIGtXMAAizeK?format=jpg&name=small](https://pbs.twimg.com/media/FSkLIGtXMAAizeK?format=jpg&name=small)

## Helper


### Helper Source code

file `MyHelper.java`

```java
import com.github.lindenb.jvarkit.tools.structvar.indexcov.SwingIndexCov;
import java.util.*;
import java.awt.Color;

public class MyHelper extends SwingIndexCov.DefaultHelper {
	private final Set<String> cases = new HashSet<>(Arrays.asList("sample1,sample2,sample3,sample4".split("[,]")));
	private final Set<Integer> cases_indexes = new HashSet<>();
	@Override
	public void initialize(final String[] header) {
		super.initialize(header);
		for(int i=3;i< header.length;i++) {
			if(cases.contains(getNormName(header[i]))) cases_indexes.add(i);
			}
		}
	private String getNormName(String sn) {
		final String[] tokens = sn.split("[_]");
		for(int i=0;i< tokens.length;i++) {
			if(tokens[i].length()==7) {
				return tokens[i];
				}
			}
		return sn;
		}
	@Override
	public String getDisplayName(final String sn) {
		String s = getNormName(sn);
		if(cases.contains(s)) return "*"+s;
		return s;
		}
	@Override
	public Optional<Colors> getColor(int[] indexes,float[] values) {
		int n_cases = 0;
		int n_other = 0;
		float treshold =  0.45f;
		for(int i=0;i< indexes.length;i++) {
			final int idx = indexes[i];
			float v  = values[i];
			boolean is_cnv = v> (1.0f + treshold) || v< (1.0f - treshold);
			if(!is_cnv) continue;
			boolean is_case = cases_indexes.contains(idx);
			if(is_case) {n_cases++;}
			else {n_other++;}
			}

		if(n_cases>0 && n_other==0) return Optional.of(Colors.RED);
		return Optional.empty();
		}
	}
```

### Compiling the helper

```make
CP=/home/lindenb/src/jvarkit-git/dist/swingindexcov.jar

myhelper.jar: MyHelper.java $(CP)
	rm -rf tmp
	mkdir -p tmp
	javac -d tmp -cp "$(CP):." $<
	jar cvf $@ -C tmp .
	rm -rf tmp
```


