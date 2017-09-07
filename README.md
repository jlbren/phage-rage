# phage-rage 1.0  

Advances in sequencing technology provide the opportunity to explore viral diversity in a variety of ecological niches. Analysis of complex viral samples presents unique bioinformatic challenges, punctuated by the paucity of available viral genomes relative to their richness in the environment. PhageRage has been developed for the comprehensive examination of viral metagenomic data. PhageRage is a UNIX-based solution for the comprehensive examination of viral metagenomic data from assembly through analysis. It was developed with agility in mind; tools, databases, and modules can easily be customized and controlled by the user without the need for extensive computational experience.

PhageRage is a locally run utility. The following README describes the process for manual installation and usage.  

PhageRage can also be quickly installed as a docker container with all dependencies included. See [here](https://github.com/thatzopoulos/PhageRage-docker) for instructions on how to build the image manually, or pull a prebuilt container from DockerHub. 
  
## Dependencies  
1. Operating System: Unix  
2. GCC Compiler: http://gcc.gnu.org/  
3. Python - Version 3.4+ 
4. KronaTools: https://github.com/marbl/Krona/wiki/KronaTools     
5. Emboss (getorf): http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html  
6. Sickle: https://github.com/ucdavis-bioinformatics/sickle *  
   1. Zlib: http://www.zlib.net/  
7. SPAdes 3.10.1: http://bioinf.spbau.ru/content/spades-download **  
8. Velvet: https://github.com/dzerbino/velvet/tree/master **  
9. Diamond: https://github.com/bbuchfink/diamond **  
10. Lambda: https://seqan.github.io/lambda/ **  
11. Megahit: https://github.com/voutcn/megahit **   
12. BLAST+ Binaries: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ **   

\* Only required when passing the --quality_control (-q) flag at runtime.   
** Only the specified assembly/alignment tool(s) selected at runtime are necesarry. 
