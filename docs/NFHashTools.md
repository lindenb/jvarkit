# NFHashTools

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

A tool to print and compare hashes .nextflow.log files


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar nfhashtools  [options] Files

Usage: nfhashtools [options] Files
  Options:
    -g, --grep
      Search for all those (sub)strings in 'hash' and 'value'. Case sensible . 
      The item matching all the regex are kept.
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -1
      hide items unique to file1
      Default: false
    -2
      hide items unique to file2
      Default: false
    -3
      hide common items in  file1 and file2
      Default: false
    -t
      look only at 'value' matching that string
    -u
      hide common fields of the unique item in  file1 and file2
      Default: false

```


## Keywords

 * json
 * nf
 * nextflow



## Creation Date

20251119

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/nextflow/nfhashtools/NFHashTools.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/nextflow/nfhashtools/NFHashTools.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/nextflow/nfhashtools/NFHashToolsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/nextflow/nfhashtools/NFHashToolsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **nfhashtools** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Motivation

this tool looks  for the log file of nextflow ( `.nextflow.log') and generate a json output.

Requirement: Nextflow run **must** be launched with `-dump-hashes json`

If there is no file, the input is read on stdin.

If there is only one file, the hashes are collected, filtered and printed as JSON.

If there are two files. Both files are read and filtered. The JSON output contains
- a list of objects uniques to the file1, 
- a list of objects uniques to the file2
- a list of common objects common to file1 and file2

If there only one object remaining in each list, an list of the common/discordant fields is also printed.


# Example


```
$ sdiff newflow.1.log newflow.2.log
Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processo   Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processo
    {                                                               {
        "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",           |         "hash": "ebc58ab2cb4848d04ec23d83f7ddf985",
        "type": "java.util.UUID",                                       "type": "java.util.UUID",
        "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"                 "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
    },                                                              },
    {                                                               {
        "hash": "f88ed0fa50638ac1efa458e70e21ebcc",                     "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
        "type": "java.lang.String",                                     "type": "java.lang.String",
        "value": "Hello"                                      |         "value": "Bonjour"
    }                                                               }
]                                                               ]
Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processo   Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processo
    {                                                               {
        "hash": "1c72b78fd00c568c524095ff174b1a6f",                     "hash": "1c72b78fd00c568c524095ff174b1a6f",
        "type": "java.util.UUID",                                       "type": "java.util.UUID",
        "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"                 "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
    },                                                              },
    {                                                               {
        "hash": "cb2d18e15fc931db472229e671a76351",                     "hash": "cb2d18e15fc931db472229e671a76351",
        "type": "java.lang.String",                                     "type": "java.lang.String",
        "value": "World"                                                "value": "World"
    }                                                               }
]                                                               ]

```

one file:

```
$ java -jar jvarkit.jar nfhashtools nextflow.1.log 2> /dev/null 
[
  [
  {
    "hash": "1c72b78fd00c568c524095ff174b1a6f",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "cb2d18e15fc931db472229e671a76351",
    "type": "java.lang.String",
    "value": "World"
  }
],
  [
  {
    "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
    "type": "java.lang.String",
    "value": "Hello"
  }
]
]
```

one file with a filter:
```
$ java -jar jvarkit.jar nfhashtools -t Hello nextflow.1.log 2> /dev/null 
[
  [
  {
    "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
    "type": "java.lang.String",
    "value": "Hello"
  }
]
]
```

with two files 

```
$ java -jar jvarkit.jar nfhashtools  nextflow.1.log nextflow.2.log 2> /dev/null 
[
  {
    "name": "L1",
    "description": "Unique to nextflow.1.log",
    "file": "nextflow.1.log",
    "items": [
      [
  {
    "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
    "type": "java.lang.String",
    "value": "Hello"
  }
]
    ]
  },
  {
    "name": "L2",
    "description": "Unique to nextflow.2.log",
    "file": "nextflow.2.log",
    "items": [
      [
  {
    "hash": "ebc58ab2cb4848d04ec23d83f7ddf985",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
    "type": "java.lang.String",
    "value": "Bonjour"
  }
]
    ]
  },
  {
    "name": "L3",
    "description": "common to nextflow.1.log and nextflow.2.log",
    "file": [
      "nextflow.1.log,
      "nextflow.2.log"
    ],
    "items": [
      [
  {
    "hash": "1c72b78fd00c568c524095ff174b1a6f",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "cb2d18e15fc931db472229e671a76351",
    "type": "java.lang.String",
    "value": "World"
  }
]
    ]
  },
  {
    "name": "uniq",
    "description": "common/discordant values shared between the unique item of L1 and the unique item of L2",
    "common": [],
    "uniq_1": [
      {
  "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
  "type": "java.lang.String",
  "value": "Hello"
},
      {
  "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",
  "type": "java.util.UUID",
  "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
}
    ],
    "uniq_2": [
      {
  "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
  "type": "java.lang.String",
  "value": "Bonjour"
},
      {
  "hash": "ebc58ab2cb4848d04ec23d83f7ddf985",
  "type": "java.util.UUID",
  "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
}
    ]
  }
]
```


