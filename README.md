# phage-rage v1.0  

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
12. [BLAST+ Binaries](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)**   

\* Only required when passing the --quality_control (-q) flag at runtime.   
\** Only the specified assembly/alignment tool(s) selected at runtime are necesarry. 


## Install/Setup
Install all necessary dependencies and make sure they are detected in the system executable search path ([PATH enviornment variable](https://help.ubuntu.com/community/EnvironmentVariables#Persistent_environment_variables)). Note that only the assemblers/mappers desired for use at runtime need to be installed. 

The required python modules can be installed quickly via pip: `pip3 install NumPy, biopython, argparse`

Dependencies and selected utilities will be checked by the pipeline at the start of each run, and an error notifying which dependency could not be located will be produced if it is unable to be found.  

To install PhageRage simply clone this git project in the desired install location. Run the pipeline 
by opening a terminal in the project directory and calling `python3 virusland.py ...`, or set up an alias/softlink as  desired. 


## Usage 
**virusland [-h] [-s | -p] [-a {spades,velvet,megahit} | -A] [-q] -m  
                 {blastp,lambda,diamond} -i INDEX [-t THREADS]  
                 [--threshold THRESHOLD] [-o OUTPUT]  
                 finput [finput ...]**  



**positional arguments:**  
  finput                
  * Input file(s). Specify 2 paired-end reads, 1 sindle-end read or 1 assembled  
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
  * Number of threads (Default=1).  

--threshold THRESHOLD  
  * Bitscore threshold for parsing mapper hits (Default=80).  

-o OUTPUT, --output OUTPUT  
  * Base output directory path. All output will be located
    here. Defaults to current working directory.  


## Examples

PhageRage can be used for the analysis of either paired end reads, single end reads, or assembled contigs. 

The following provide sample runs for each case:

**Paired-end reads:** 

```virusland.py /path/to/read1.fa /path/to/read2.fa -pqa spades -m diamond -i /path/to/GBK_files -t 12 -o /my/output/dir ``` 
  
or 

```virusland.py /path/to/read1.fa /path/to/read2.fa --paired_end_reads --quality_control --assembler spades --mapper diamond --index /path/to/GBK_files --threads --output /my/output/dir```  


The above command will run the pipeline on 2 paired end reads using 12 threads, perform quality control timming with Sickle, assembly with Spades, and map reads with diamonds. 

**Single-end reads:**  

```virusland.py /path/to/myread.fa -sa velvet -m blast -i /path/to/GBK_files -o /my/output/dir ```

or

```virusland.py /path/to/myread.fa --single_reads --assembler velvet --mapper blast --index /path/to/GBK_files ---output /my/output/dir``` 

The above command will run the pipeline on single-end reads using a default of 1 thread, perform assembly with Velvet, and mapping with BLAST. 

**Assembled contigs:**

```virusland /path/to/my_contigs.fa -Am lambda -i /path/to/GBK_files --threshold .8 -o /my/output/dir ```

or

```virusland /path/to/my_contigs.fa --assembled_contigs --mapper lambda --index /path/to/GBK_files --threshold 60 --output /my/output/dir ``` 

The above command accepts a single contigs file and performs mapping using lambda with a bitscore threshold of 60.


## Output

All output will be located in the following respective subdirectories within the base output directory specified by the -o/--output flag. Note that if any of the created subdirectories already exist in the specified output location an error will be raised in order avoid over writting previous runs. 

**Stats**

Subdirectory containing the final output and statistics produced by the pipeline.
 * coverage.csv: File providing stats on total percentage of each organism's genome was mapped to hits.
 * hits_by_protein.csv: File providing counts of how many hits each individual protein received.
 * krona_stats.csv: File produced that can be turned into a Krona HTML graph by the ktImportText utility. 
 * krona_graph.html: Krona hierarchial data graph. Can be opened by any modern web browser. 
 * hitviz_stats.csv: File produced for use with the [hitviz](<hitviz_repo>) visualization utility. 
 
**Trimmed**
 
 Subdirectory containing trimmed produced by Sickle if the quality control flag is passed at runtime.  
  * trimmed-myread.fastq .. : Trimmed read file(s) produced by Sickle. 
  * singletons.fastq: Unpaired reads produced by Sickle when processing paired-end reads. 
  
  
**Assembled**

Subdirectory containing the output of the selected assembly utility.  
 * contigs.fasta: Contigs file produced by the assembler.
  
**Mapped**
  
  Subdirectory containing mapper output and related files.  
  * hits.csv: Hits file produced by selected mapper utility. 
  * index.faa: FAA file produced by parsing the GBK files provided by the -i/--index flag at runtime used to build the mapper index. 
  * predicted_orfs.faa: Proteins predicted from the provided or assembled contigs produced by getorf that are used for mapping against the index. 
  
  **Logs**
  
  * Directory containing pipeline log file which records user supplied parameters, runtime per module, and additional supplementary information. 

