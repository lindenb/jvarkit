# UKBiobankSelectSamples

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Select samples from ukbiobank


## Usage

```
Usage: ukbiobanksamples [options] Files
  Options:
    --column
      column pattern. Look for columns starts with 'x'
      Default: f.41270.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hide-phenotypes
      do not print phenotypes
      Default: false
    --inverse
      inverse selection of sample when using option -A. get samples NOT having 
      the phenotypes
      Default: false
  * --ontology, --coding
      coding/ontology file.
    -o, --output
      Output file. Optional . Default: stdout
  * --tab
      *.tab file.
    --user-coding, -A
      limit to those coding. multiple separated by commas
      Default: []
    --version
      print version and exit

```


## Keywords

 * ukbiobank


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew ukbiobanksamples
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210705

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ukbiobank/UKBiobankSelectSamples.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ukbiobank/UKBiobankSelectSamples.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **ukbiobanksamples** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/ukbiobanksamples.jar  --ontology /path/to/coding19.tsv --tab  /path/to/ukb46178.tab --column "f.41270." -A 'I34'

[INFO][UKBiobankSelectSamples]using I349{I34.9 Nonrheumatic mitral valve disorder, unspecified} I34{I34 Nonrheumatic mitral valve disorders} I340{I34.0 Mitral (valve) insufficiency} I341{I34.1 Mitral (valve) prolapse} I342{I34.2 Nonrheumatic mitral (valve) stenosis} I348{I34.8 Other nonrheumatic mitral valve disorders}

1000767	R55 [R55 Syncope and collapse ]	W190 [W19.0 Home ]	Z961 [Z96.1 Presence of intraocular lens ]	K529 [K52.9 Non-infective gastro-enteritis and colitis, unspecified ]	K638 [K63.8 Other specified diseases of intestine ]	H041 [H04.1 Other disorders of lachrymal gland ]	I739 [I73.9 Peripheral vascular disease, unspecified ]	H269 [H26.9 Cataract, unspecified ]	R104 [R10.4 Other and unspecified abdominal pain ]	N946 [N94.6 Dysmenorrhoea, unspecified ]	N905 [N90.5 Atrophy of vulva ]	E119 [E11.9 Without complications ]	N908 [N90.8 Other specified noninflammatory disorders of vulva and perineum ]	N920 [N92.0 Excessive and frequent menstruation with regular cycle ]	I341 [I34.1 Mitral (valve) prolapse ] (...)
1002786	I10 [I10 Essential (primary) hypertension ]	R11 [R11 Nausea and vomiting ]	R12 [R12 Heartburn ]	K20 [K20 Oesophagitis ]	Z864 [Z86.4 Personal history of psychoactive substance abuse ]	M159 [M15.9 Polyarthrosis, unspecified ]	Z423 [Z42.3 Follow-up care involving plastic surgery of upper extremity ]	W010 [W01.0 Home ]	I120 [I12.0 Hypertensive renal disease with renal failure ]	J459 [J45.9 Asthma, unspecified ]	K299 [K29.9 Gastroduodenitis, unspecified ]	D509 [D50.9 Iron deficiency anaemia, unspecified ]	D125 [D12.5 Sigmoid colon ]	E835 [E83.5 Disorders of calcium metabolism ]	K298 [K29.8 Duodenitis ]	M170 [M17.0 Primary gonarthrosis, bilateral ]	I340 [I34.0 Mitral (valve) insufficiency ]...
1003404	R02 [R02 Gangrene, not elsewhere classified ]	E86 [E86 Volume depletion ]	J40 [J40 Bronchitis, not specified as acute or chronic ]	Z602 [Z60.2 Living alone ]I420 [I42.0 Dilated cardiomyopathy ]	I071 [I07.1 Tricuspid insufficiency ]	I340 [I34.0 Mitral (valve) insufficiency ]
1003703	M545 [M54.5 Low back pain ]	M754 [M75.4 Impingement syndrome of shoulder ]	I340 [I34.0 Mitral (valve) insufficiency ](...)
```

