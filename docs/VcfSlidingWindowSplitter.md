# VcfSlidingWindowSplitter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split VCF by sliding window


## Usage

```
Usage: vcfwindowsplitter [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --manifest
      Manifest Bed file output containing chrom/start/end of each gene
    -M, --max-variant
      Maximum number of variants required to write a vcf. don't write if 
      num(variant) > 'x' . '<=0' is ignore
      Default: -1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -n, --min-variant
      Minimum number of variants required to write a vcf. don't write if 
      num(variant) < 'x'
      Default: 1
  * -o, --output
      An existing directory or a filename ending with the '.zip' suffix.
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit
    -s, -S, --window-shift
      Sliding window shift. A distance specified as a positive integer.Commas 
      are removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 500000
    -w, -W, --window-size
      Sliding window size. A distance specified as a positive integer.Commas 
      are removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 1000000

```


## Keywords

 * vcf
 * sliding
 * window


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfwindowsplitter
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190619

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfSlidingWindowSplitter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfSlidingWindowSplitter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfSlidingWindowSplitterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfSlidingWindowSplitterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfwindowsplitter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example

```
$ java -jar dist/vcfwindowsplitter.jar -n 2 -w 1000 -s 500 -o jeter.zip -m jeter.manifest src/test/resources/rotavirus_rf.vcf.gz 
[INFO][VcfSlidingWindowSplitter]. Completed. N=45. That took:0 second

$ head jeter.manifest  | column -t
#chrom  start  end   window          path                                                     Count_Variants
RF02    250    877   RF02:0-1000     35/a02abf73216a95f62458f41dfbbf79/RF02_0_1000.vcf.gz     3
RF02    577    877   RF02:500-1500   b5/ed3c205185d506e382c437cc010dee/RF02_500_1500.vcf.gz   2
RF02    1725   1965  RF02:1000-2000  41/026109d82e6f553dd0b6244ceffe5a/RF02_1000_2000.vcf.gz  2
RF02    1725   1965  RF02:1500-2500  2f/82d6416c74c8b2427fb71a696bf491/RF02_1500_2500.vcf.gz  2
RF03    1220   1242  RF03:500-1500   3b/11b08845e1d70e4b4c6ab3a5a279ee/RF03_500_1500.vcf.gz   2
RF03    1220   1708  RF03:1000-2000  dc/e6069ebdb06ef2cfd85e96331d29d2/RF03_1000_2000.vcf.gz  4
RF03    1687   2315  RF03:1500-2500  af/b0daefa93f1423494ed6da29ab49c7/RF03_1500_2500.vcf.gz  5
RF03    2149   2573  RF03:2000-3000  f7/ffe3fb2204dd071e7d40d8e0e1151e/RF03_2000_3000.vcf.gz  4
RF04    886    991   RF04:0-1000     5b/9e9699f08473c9fa7cd550247c36ff/RF04_0_1000.vcf.gz     2

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     1359  2019-06-19 14:48   35/a02abf73216a95f62458f41dfbbf79/RF02_0_1000.vcf.gz
     1283  2019-06-19 14:48   b5/ed3c205185d506e382c437cc010dee/RF02_500_1500.vcf.gz
     1301  2019-06-19 14:48   41/026109d82e6f553dd0b6244ceffe5a/RF02_1000_2000.vcf.gz
     1301  2019-06-19 14:48   2f/82d6416c74c8b2427fb71a696bf491/RF02_1500_2500.vcf.gz
     1281  2019-06-19 14:48   3b/11b08845e1d70e4b4c6ab3a5a279ee/RF03_500_1500.vcf.gz
     1392  2019-06-19 14:48   dc/e6069ebdb06ef2cfd85e96331d29d2/RF03_1000_2000.vcf.gz
     1462  2019-06-19 14:48   af/b0daefa93f1423494ed6da29ab49c7/RF03_1500_2500.vcf.gz
     1418  2019-06-19 14:48   f7/ffe3fb2204dd071e7d40d8e0e1151e/RF03_2000_3000.vcf.gz
     1265  2019-06-19 14:48   5b/9e9699f08473c9fa7cd550247c36ff/RF04_0_1000.vcf.gz
     1405  2019-06-19 14:48   4c/75302cbe2cd7f16759c5e864031e03/RF04_500_1500.vcf.gz
     1482  2019-06-19 14:48   d4/7c942d8e574e6bdfccafbeb767675e/RF04_1000_2000.vcf.gz
     1356  2019-06-19 14:48   7f/af3063f079076117e4186fab17768b/RF04_1500_2500.vcf.gz
     1414  2019-06-19 14:48   e6/73552f774e54d0c0bd8bfc31192335/RF05_0_1000.vcf.gz
     1410  2019-06-19 14:48   95/3b7f4b269e5dc25550ed7610c4be52/RF05_500_1500.vcf.gz
     1274  2019-06-19 14:48   1e/a8b7283d053ccf22f44301504d522f/RF05_1000_2000.vcf.gz
     1427  2019-06-19 14:48   ea/ef4a24f964402115bfbf6b9f25253e/RF06_0_1000.vcf.gz
     1504  2019-06-19 14:48   d5/487160765c2bc2a1d286e534c41d06/RF06_500_1500.vcf.gz
     1406  2019-06-19 14:48   15/c37e21e6933c610a8bf53d833fcea9/RF07_0_1000.vcf.gz
     1270  2019-06-19 14:48   ff/af5bab3250e5ea8106efc4ca424e4f/RF07_500_1500.vcf.gz
     1269  2019-06-19 14:48   4b/9c386cb83786e3a61bbc54b949d5f9/RF08_0_1000.vcf.gz
     1270  2019-06-19 14:48   bb/fc42ed842a1d5f76bd19986655c7f7/RF08_500_1500.vcf.gz
     1333  2019-06-19 14:48   07/9587be0df182f9e938d412ffd9ffa6/RF09_0_1000.vcf.gz
     1341  2019-06-19 14:48   b1/5a51dc37c5b128a2edd36a68fcb9cb/RF10_0_1000.vcf.gz
---------                     -------
    31223                     23 files

```

