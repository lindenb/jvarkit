JVARKIT
=======

Author      : Pierre Lindenbaum Phd. Institut du Thorax. Nantes. France.
Version     : 8f0d6d899
Compilation : 20260131015012
Github      : https://github.com/lindenb/jvarkit
Issues      : https://github.com/lindenb/jvarkit/issues

## Usage

```
  java -jar jvarkit.jar [options]
```
or
```
  java -jar jvarkit.jar <command name> (other arguments)
```

## Options

 + --help show this screen
 + --help-all show all commands, including the private ones.
 + --version print version

## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download a pre-compiled executable jar

A pre-compiled java executable of jvarkit.jar is available at [https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH](https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH)
### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew jvarkit
```

The java jar file will be installed in the `dist` directory.

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **jvarkit** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

