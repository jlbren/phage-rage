# phage-rage 1.0  

Advances in sequencing technology provide the opportunity to explore viral diversity in a variety of ecological niches. Analysis of complex viral samples presents unique bioinformatic challenges, punctuated by the paucity of available viral genomes relative to their richness in the environment. PhageRage has been developed for the comprehensive examination of viral metagenomic data. PhageRage is a UNIX-based solution for the comprehensive examination of viral metagenomic data from assembly through analysis. It was developed with agility in mind; tools, databases, and modules can easily be customized and controlled by the user without the need for extensive computational experience.

PhageRage is a locally run utility. The following README describes the process for manual installation and usage.  

PhageRage can also be quickly installed as a docker container with all dependencies included. See [here](https://github.com/thatzopoulos/PhageRage-docker) for instructions on how to build the image manually, or pull a prebuilt container from DockerHub. 
  
## Dependencies  
1. Operating System: Unix  
2. [GCC Compiler](http://gcc.gnu.org/)  
3. Python - Version 3.4+   
    * [Numpy](https://www.scipy.org/scipylib/download.html)
    * [Biopython](http://biopython.org/wiki/Download)
    * ArgParse
4. [KronaTool](https://github.com/marbl/Krona/wiki/KronaTools)    
5. [Emboss (getorf)](http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html)  
6. [Sickle](https://github.com/ucdavis-bioinformatics/sickle) \*  
    * [Zlib](http://www.zlib.net/)  
7. [SPAdes v3.10.1](http://bioinf.spbau.ru/content/spades-download) \**  
8. [Velvet](https://github.com/dzerbino/velvet/tree/master) \**  
9. [Diamond](https://github.com/bbuchfink/diamond) \**  
10. [Lambda](https://seqan.github.io/lambda/) \**  
11. [Megahit](https://github.com/voutcn/megahit) \**   
12. [BLAST+ Binaries](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) \**   

\* Only required when passing the --quality_control (-q) flag at runtime.   
\** Only the specified assembly/alignment tool(s) selected at runtime are necesarry. 


## Install/Setup
Install all necessary dependencies and make sure they are detected in the system executable search path ([PATH enviornment variable](https://help.ubuntu.com/community/EnvironmentVariables#Persistent_environment_variables)). Note that only the assemblers/mappers desired for use at runtime need to be installed. 

The required python modules can be installed quickly via pip: `pip3 install NumPy, biopython, argparse`

Dependencies and selected utilities will be checked by the pipeline at the start of each run, and an error notifying which dependency could not be located will be produced if it is unable to be found.  


## Usage 
**virusland [-h] [-s | -p] [-a {spades,velvet,megahit} | -A] [-q] -m  
                 {blastp,lambda,diamond} -i INDEX [-t THREADS]  
                 [--threshold THRESHOLD] [-o OUTPUT]  
                 finput [finput ...]**  



**positional arguments:**  
  finput                
  * Input file(s). Specify 2 PE, 1 SE or 1 assembled  
    contigs file(s).  
  
**optional arguments:**  
-h, --help            
  * Show this help message and exit  

-s, --single_reads    
  * Single reads flag. Requires specifying 1 input read file.  

-p, --paired_end_reads  
  * Paired-end reads flag. Requires specifying 2 input read files.  

-a {spades,velvet,megahit}, --assembler {spades,velvet,megahit}  
  * Selected assembly utility.  

-A, --assembled_contigs  
  * Assembled contigs flag. Start analysis from an
    existing assembly. Requires specifying 1 input contigs
    file.  

-q, --quality_control  
  * Read quality control flag. Trim input read(s) with sickle.  

-m {blastp,lambda,diamond}, --mapper {blastp,lambda,diamond}  
  * Selected mapping utility.  

-i INDEX, --index INDEX  
  * Path to directory containing GBK files to be used when
    building mapper index. Files may either be within the
    given directory and/or one level of subdirectory.

-t THREADS, --threads THREADS  
  * Number of threads.  

--threshold THRESHOLD  
  * Bitscore threshold for parsing mapper hits.  

-o OUTPUT, --output OUTPUT  
  * Base output directory path. All output will be located
    here. Must be relative to current working directory.  
