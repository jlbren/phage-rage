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
4. [KronaTools](https://github.com/marbl/Krona/wiki/KronaTools)    
5. [Emboss (getorf)](http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html)  
6. [Sickle](https://github.com/ucdavis-bioinformatics/sickle) \*  
    * [Zlib](http://www.zlib.net/)  
7. [SPAdes v3.10.1](http://bioinf.spbau.ru/content/spades-download) \**  
8. [Velvet](https://github.com/dzerbino/velvet/tree/master) \**  
9. [Diamond](https://github.com/bbuchfink/diamond) \**  
10. [Lambda](https://seqan.github.io/lambda/) \**  
11. [Megahit](https://github.com/voutcn/megahit) \**   
12. [BLAST+ Binaries](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) **   

\* Only required when passing the --quality_control (-q) flag at runtime.   
\** Only the specified assembly/homology detection tool(s) selected at runtime are necessary. 


## Install/Setup
Install all necessary dependencies and make sure they are detected in the system executable search path ([PATH environment variable](https://help.ubuntu.com/community/EnvironmentVariables#Persistent_environment_variables)). Note that only the assembly/homology detection tool desired for use at runtime need to be installed. 

The required python modules can be installed quickly via pip: `pip3 install NumPy, biopython, argparse`

Dependencies and selected utilities will be checked by the pipeline at the start of each run, and an error notifying which dependency could not be located will be produced if it is unable to be found.  

To install PhageRage simply clone this git project in the desired install location. Run the pipeline 
by opening a terminal in the project directory and calling `python3 phage-rage.py ...`, or set up an alias/softlink as  desired. 


## Usage 
**phage-rage [-h] [-s | -p] [-a {spades,velvet,megahit} | -A] [-q] -m  
                 {blastp,lambda,diamond} -i INDEX [-t THREADS]  
                 [--threshold THRESHOLD] [-o OUTPUT]  
                 finput [finput ...]**  



**positional arguments:**  
  finput                
  * Input file(s). Specify 2 paired-end reads, 1 single-end read, or 1 assembled  
    contigs file(s).  
  
**optional arguments:**  
-h, --help            
  * Show this help message and exit  

-s, --single_reads    
  * Single-end reads flag. Requires specifying 1 input read file.  

-p, --paired_end_reads  
  * Paired-end reads flag. Requires specifying 2 input read files.  

-a {spades,velvet,megahit}, --assembler {spades,velvet,megahit}  
  * Selected assembly utility.  

-A, --assembled_contigs  
  * Assembled contigs flag. Start analysis from an
    existing assembly. Requires specifying 1 input contig
    file.  

-q, --quality_control  
  * Read quality control flag. Trim input read(s) with Sickle.  

-m {blastp,lambda,diamond}, --mapper {blastp,lambda,diamond}  
  * Selected homology detection tool.  

-i INDEX, --index INDEX  
  * User supplied path to directory containing GBK files to be used by
    the homology detection tool. Files may either be within the
    specified directory and/or one level of subdirectory.

-t THREADS, --threads THREADS  
  * Number of threads. (Default=1) 

--threshold THRESHOLD  
  * Bitscore threshold for parsing mapper hits. (Default=80) 

-o OUTPUT, --output OUTPUT  
  * Base output directory path. All output will be located
    here. (Default=current working directory)  


## Examples

PhageRage can be used for the analysis of either paired end reads, single-end reads, or assembled contigs. 

The following provide sample runs for each case:

**Paired-end reads:** 

```phage-rage.py /path/to/read1.fa /path/to/read2.fa -pqa spades -m diamond -i /path/to/GBK_files -t 12 -o /my/output/dir ``` 
  
or 

```phage-rage.py /path/to/read1.fa /path/to/read2.fa --paired_end_reads --quality_control --assembler spades --mapper diamond --index /path/to/GBK_files --threads --output /my/output/dir```  


The above command will run the pipeline on 2 paired end reads using 12 threads, perform quality control trimming with Sickle, assembly with Spades, and homology detection with Diamond. 

**Single-end reads:**  

```phage-rage.py /path/to/myread.fa -sa velvet -m blast -i /path/to/GBK_files -o /my/output/dir ```

or

```phage-rage.py /path/to/myread.fa --single_reads --assembler velvet --mapper blast --index /path/to/GBK_files ---output /my/output/dir``` 

The above command will run the pipeline on single-end reads using a default of 1 thread, perform assembly with Velvet, and homology detection with blastp. 

**Assembled contigs:**

```phage-rage.py /path/to/my_contigs.fa -Am lambda -i /path/to/GBK_files --threshold .8 -o /my/output/dir ```

or

```phage-rage.py /path/to/my_contigs.fa --assembled_contigs --mapper lambda --index /path/to/GBK_files --threshold 60 --output /my/output/dir ``` 

The above command accepts a contigs file and performs homology detection using Lambda with a bitscore threshold of 60.


## Output

All output will be located in the following respective subdirectories within the base output directory specified by the -o/--output flag. Note: if any of the created subdirectories already exist in the specified output location, an error will be raised in order avoid overwriting previous runs.

**Stats**

Subdirectory containing the final output and statistics produced by the pipeline.
 * coverage.csv: File providing stats on total percentage of each sequence in the index set identified as homologous to contigs.
 * hits_by_protein.csv: File providing counts of hits for each individual coding sequence in each sequence in the index set.
 * krona_stats.csv: File produced that can be turned into a Krona HTML graph by the ktImportText utility. 
 * krona_graph.html: Krona hierarchical data graph. Can be opened by any modern web browser. 
 * hitviz_stats.csv: File produced for use with the [hitviz](<hitviz_repo>) visualization utility. 
 
**Trimmed**
 
 Subdirectory containing trimmed produced by Sickle if the quality control flag is passed at runtime.  
  * trimmed-myread.fastq .. : Trimmed read file(s) produced by Sickle. 
  * singletons.fastq: Unpaired reads produced by Sickle when processing paired-end reads. 
  
  
**Assembled**

Subdirectory containing the output of the selected assembler.  
 * contigs.fasta: Contigs file produced by the assembler.
  
**Mapped**
  
  Subdirectory containing output from the selected homology detection tool and related files.  
  * hits.csv: Hits file produced by selected homology detection tool. 
  * index.faa: FAA file produced at runtime by parsing the GBK files provided by the -i/--index flag. 
  * predicted_orfs.faa: Proteins predicted from the provided or assembled contigs produced by getorf. These sequences are used for querying against the index. 
  
  **Logs**
  
  * Directory containing PhageRage log file which records user supplied parameters, runtime per module, and additional supplementary information. 
  
  ## Authors
  
  * Jonathon Brenner
  * Thomas Hatzopoulos
  * Zach Romer 
  * Catherine Putonti 
  
  ### Contact: putonti.lab@gmail.com.

