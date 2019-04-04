# Optimizer

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Genetic-Programming-like parameters optimizer


## Usage

```
Usage: optimizer [options] Files
  Options:
    -A, --all
      Run all possible combinations
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -seed, --random
      Random seed. -1 == current time
      Default: -1
  * -c, --code, --source
      User's code.
    --version
      print version and exit

```


## Keywords

 * genetic-programming


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew optimizer
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/optimizer/Optimizer.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/optimizer/Optimizer.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **optimizer** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


This is a Genetic-Programming-like parameters optimizer.

use must provide the scafold of a java code that will override `Solution` and a JSON file describing the
parameters.

## An example of json config

```json

{
"params":[
	{
	"name":"param1",
	"type":"int",
	"min": 1,
	"max":10,
	"shift":1
	},
	{
	"name":"param2",
	"type":"double",
	"min": 0.01,
	"max": 0.1,
	"shift":0.01
	}
  ]
}
```


## The base class Solution

```java
public static abstract class Solution implements Comparable<Solution>
	{
	protected long generation =  -1L;
	protected final Map<String,Object> params;
	public Solution(final Map<String,Object> params)
		{
		this.params = Collections.unmodifiableMap(params);
		}
	// eval the result. Must be implemented by the user
	public abstract int execute() throws Exception;
	// delete any file associated to this solution 
	public void delete() {}
	//mate a Solution with another, returns null if mating is not possible
	public Solution mate(final Solution another) {
		return null;
	}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj ==null || !(obj instanceof Solution)) return false;
		final Solution other=Solution.class.cast(obj);
		return this.params.equals(other.params);
		}
	
	@Override
	public int hashCode() {
		return this.params.hashCode();
		}
	}
```

## An example of custom class extending `Solution`

`__BASE__` will be replaced by the base class name `Solution`.
`__CLASS__` will be replaced by the current generated class name.

The user's code will be inserted in the following template:

```
 1  import java.util.*;
 2  import java.io.*;
 3  import java.util.stream.*;
 4  import java.util.function.*;
 5  import htsjdk.samtools.util.*;
 6  import htsjdk.variant.variantcontext.*;
 7  import htsjdk.variant.vcf.*;
 8  import javax.annotation.processing.Generated;
 9  @Generated(value="Optimizer",date="2017-07-10T11:20:07+0200")
10  public class __CLASS__ extends  __BASE__ {
11  public __CLASS__(final Map<String,Object> params) {
12  super(params);
13  }
14      // user's code starts here 
(...)  
93     // user's code ends here 
94  }

```

in __CLASS__ User must implement:

* 'compareTo' to compare two solutions
* 'execute' to compute the result with the current params. Returns '0' on success.
* 'delete' remove resources associated to this Solution.

