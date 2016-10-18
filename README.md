# SNPprimers

# Requirements

If you plan on running installing the required modules using setup.py, following an `sudo apt-get update` the following packages need to be installed (these are the packages that need to be installed in an Ubuntu:16.04 docker container):

* git
* python-dev
* python-pip

`sudo apt-get install -y build-essential git python-dev python-pip` 

# Installation

Clone the repository (--recursive will clone the necessary submodules):

`git clone https://github.com/adamkoziol/SNPprimers.git --recursive`

Install python dependencies:

	
```
cd SNPprimers/
python setup.py install
```

# Usage
## snpprimers.py

Extracts nucleotide sequence data surrounding SNP coordinates in .fasta assembly files. Create primers to amplify extracted regions using Primer3

### Example command

`python snpprimers.py -c coords.txt /path/to/desired/output/folder`

Required arguments:

* path to the folder to be used to store results
* -c: File with .txt extension containing one contig name and SNP coordinate pair separated by a delimiter (see below) per line. If a path is not specified in the file name, then the program will assume that the file is in the path supplied above

Optional arguments:

* -d: The delimiter used to separate contig name and SNP coordinates. Popular options are "space", "tab", and "comma". Default is space. Be aware that a delimiter, 'such as "-" will break the program if there are hyphens in your contig names
* -a: The size of the amplicon to extract from the assembly file. The default is 300; 150 bp upstream and 150 bp downstream of the SNP
* -s: The location of the .fasta assembly file. If not supplied, the program will assume that this file is in path/sequences


See usage below:

```
usage: snpprimers.py [-h] -c COORDINATEFILE [-d DELIMITER] [-a AMPLICONSIZE]
                     [-s SEQUENCEPATH]
                     path

Extracts nucleotide sequence data surrounding SNP coordinates in .fasta
assembly files.Create primers to amplify extracted regions using Primer3

positional arguments:
  path                  Specify path of folder in which the analyses are to be
                        performed

optional arguments:
  -h, --help            show this help message and exit
  -c COORDINATEFILE, --coordinateFile COORDINATEFILE
                        File with .txt extension containing one contig name
                        and SNP coordinate pair separated by a delimiter (see
                        below) per line. If a path is not specified in the
                        file name, then the program will assume that the file
                        is in the path supplied above.
  -d DELIMITER, --delimiter DELIMITER
                        The delimiter used to separate contig name and SNP
                        coordinates. Popular options are "space", "tab", and
                        "comma". Default is space. Be aware that a delimiter,
                        such as "-" will break the program if there are
                        hyphens in your contig names
  -a AMPLICONSIZE, --ampliconSize AMPLICONSIZE
                        The size of the amplicon to extract from the assembly
                        file. The default is 300; 150 bp upstream and 150 bp
                        downstream of the SNP. I suppose the amplicon will
                        actually be 301 bp, since the SNP position needs to be
                        included, but round numbers are easier
  -s SEQUENCEPATH, --sequencepath SEQUENCEPATH
                        The location of the .fasta assembly file. If not
                        supplied, the program will assume that this file is in
                        path/sequences
```
