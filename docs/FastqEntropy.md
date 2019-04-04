# FastqEntropy

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compute the Entropy of a Fastq file (distribution of the length(gzipped(sequence))


## Usage

```
Usage: fastqentropy [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * fastq


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew fastqentropy
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqEntropy.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/FastqEntropy.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastqentropy** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$  java -jar dist/fastqentropy.jar file.fastq.gz
```
output (length/count):
```
36	1
37	1
38	2
39	2
128	3
41	7
40	9
127	10
42	18
43	20
44	34
126	35
45	71
46	114
125	114
47	198
124	254
48	272
49	475
123	552
50	579
51	787
52	1028
122	1085
53	1281
54	1391
55	1624
121	1846
56	1847
57	1907
58	2182
60	2234
59	2245
62	2387
61	2392
63	2460
64	2478
66	2613
65	2655
67	2740
68	2753
78	2881
77	2919
79	2921
76	2928
75	2954
80	2958
69	2969
70	3014
81	3014
71	3030
72	3056
82	3067
74	3073
73	3075
120	3157
83	3175
84	3381
85	3433
86	3627
87	3846
88	4098
89	4417
90	4535
119	4951
91	5096
92	5544
93	6045
94	6475
95	6913
118	7089
96	7211
97	7697
98	7877
99	8118
100	8194
101	8235
102	8280
103	8496
104	8740
105	9144
117	9335
106	9571
107	10767
116	11688
108	12174
109	13664
115	14461
110	15637
114	16983
111	17347
113	18332
112	18698
```
