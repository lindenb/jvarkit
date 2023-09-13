JVARKIT
=======

Java utilities for Bioinformatics

[![Build Status](https://travis-ci.org/lindenb/jvarkit.svg)](https://travis-ci.org/lindenb/jvarkit)

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

## Documentation

Documentation is available at: [https://jvarkit.readthedocs.io/](https://jvarkit.readthedocs.io/)

## Download

A pre-compiled jar is available at [https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH](https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH) . 

## Compilation

Since 2023, most tools (but not all) are now packaged into one application `jvarkit.jar`. Tools that were executed like `java -jar toolname.jar` are now executed as `java -jar jvarkit.jar toolname`. The documentation is not always up to date on this point.

See the documentation at [https://jvarkit.readthedocs.io/](https://jvarkit.readthedocs.io/).

## Conda / Bioconda

A conda package created by [DrYak](https://github.com/DrYak)  is available at [https://bioconda.github.io/recipes/jvarkit/README.html](https://bioconda.github.io/recipes/jvarkit/README.html)


## Containers

jvarkit is available as a **Docker** container at [https://hub.docker.com/r/lindenb/jvarkit](https://hub.docker.com/r/lindenb/jvarkit) . The jar is compiled under `/opt/jvarkit/dist/jvarkit.jar` so a command should be  `docker run java -jar /opt/jvarkit/dist/jvarkit.jar` . Nevertheless GUI/Swing applications don't work. 

## Author

Pierre Lindenbaum PhD

[@yokofakun](https://twitter.com/yokofakun)


