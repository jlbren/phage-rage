#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>

using namespace std;

bool pauda, lambda, other;

class genome
{
    public:
    int num_proteins;       //number of proteins in this genome
    char ** proteins;       //1D array of the protein accession numbers (names)
    int * counts;           //number of times the protein is hit (indices of counts correspond with the indices of proteins)
    char * acc_genome;      //accession # of the genome
    char * spp_information; //organism's name

    genome();
    bool initialize(int i, char * accession, char * spp);
    bool initialize_counts();
    bool add_protein(int i, char * val);
    int get_number_of_hits_to_genome();
    int get_number_of_proteins_with_hits();
    bool get_number_of_hits_by_protein(char **& protein_id, int *& vals);
    bool show_number_of_hits_by_protein();
    float get_percentage_of_proteins_hit();
    bool increment_hits_counter(int index);
    bool show_hits();
    ~genome();
};

genome::genome()
{
    num_proteins=0;
}

bool genome::initialize(int i, char * accession, char * spp)
{
    num_proteins=i;
    acc_genome=new char [strlen(accession)+1];
    sprintf(acc_genome,"%s",accession);
    spp_information=new char [strlen(spp)+1];
    sprintf(spp_information,"%s",spp);

    proteins=new char * [num_proteins];
    counts=new int [num_proteins];
    initialize_counts();
}

bool genome::initialize_counts()
{
    int i;
    for(i=0;i<num_proteins;i++) counts[i]=0;
    return true;
}

bool genome::add_protein(int i, char * val)
{
    if(i<num_proteins)
    {
        proteins[i]=new char [strlen(val)+1];
        sprintf(proteins[i],"%s",val);
        return true;
    }
    else return false;
}

int genome::get_number_of_hits_to_genome()
{
    int i;
    int total=0;
    for(i=0;i<num_proteins;i++) total+=counts[i];
    return total;
}


int genome::get_number_of_proteins_with_hits()
{
    int i;
    int counter=0;
    for(i=0;i<num_proteins;i++)
    {
        if(counts[i]>0) counter++;
    }
    return counter;
}

bool genome::get_number_of_hits_by_protein(char **& protein_id, int *& vals)
{
    if(vals!=NULL) delete vals;
    vals=new int [num_proteins];
    if(protein_id!=NULL) delete protein_id;
    protein_id=new char * [num_proteins];
    int i;
    for(i=0;i<num_proteins;i++)
    {
        vals[i]=counts[i];
        protein_id[i]=new char [strlen(proteins[i])+1];
        sprintf(protein_id[i],"%",proteins[i]);
    }
    return true;
}

bool genome::show_hits()
{
    int i;
    for(i=0;i<num_proteins;i++)
    {
        if(counts[i]>0)
        {
            cout << proteins[i] << "\t" << counts[i] << endl;
        }
    }
}

bool genome::show_number_of_hits_by_protein()
{
    int i;
    for(i=0;i<num_proteins;i++)
    {
        cout << proteins[i] << "\t" << counts[i] << endl;
    }
    return true;
}


float genome::get_percentage_of_proteins_hit()
{
    int i;
    int n=0;
    for(i=0;i<num_proteins;i++)
    {
        if(counts[i]>0) n++;
    }
    return ((float)n/(float)num_proteins)*100.00;
}

bool genome::increment_hits_counter(int index)
{
    if(index<num_proteins) {counts[index]++; return true;}
    else return false;
}

genome::~genome()
{
    int i;
    for(i=0;i<num_proteins;i++) delete proteins[i];
    delete proteins;
    delete counts;
}

class super_file
{
    public:
    int num_genomes;        //number of genomes in faa collection
    genome * g;             //array of type genome

    super_file(int n, char * file);
    bool read_BLAST(char * file);
    bool read_BLAST(char * file, float evalue_threshold);
    bool is_in(char * acc, int & genome_num, int & gene_num);
    bool write_out_summary_statistics(char * output_file, char * output_root);
    bool write_out_forHitViz(char * output_file, char * output_root);
    bool write_out_forKrona(char * output_file, char * output_root);
    bool reset_counts();
    ~super_file();

    private:
    bool populate_genome(char * file);
    bool parse_hit_special_case(char *& hit, char *& prot_acc_num, char *& prot_function, char *& spp);
    bool parse_hit(char *& hit, char *& prot_acc_num, char *& prot_function, char *& spp);
    bool parse_hit(char *& hit, char *& prot_acc_num);
};

super_file::super_file(int n, char * file)
{
    num_genomes=n;
    g=new genome [num_genomes];
    populate_genome(file);
    cout << "super_file set up\n";
    cout << "number of genomes " << num_genomes << endl;
}

bool super_file::populate_genome(char * file)
{   cout << "Populating genome" << endl;
    int i;
    char * text=new char [100];
    char * spp=new char [100];
    char * prior_g=new char [100];
    int num_genes=0;
    int genome_index=0;

    ifstream in;

    //set up data structure
    in.open(file);
    cout << file << endl;
    if(!in.is_open())
    {
        cout << "Cannot open genome file.\n";
	    cout << file << endl;
        delete text;
        delete prior_g;
        delete spp;
        in.clear();
        in.close();
        return false;
    }
    in.getline(prior_g,100,'\t'); num_genes++;
    in.getline(text,100,'\t');
    in.getline(spp,100,'\n');
    while(in.peek()!=EOF)
    {
        in.getline(text,100,'\t');
        if(strcmp(text,prior_g)==0)
            num_genes++;
        else
        {
            g[genome_index].initialize(num_genes,prior_g,spp);
            num_genes=1;
            genome_index++;
            sprintf(prior_g,"%s",text);
        }
        in.getline(text,100,'\t');
        in.getline(spp,100,'\n');
    }
    g[genome_index].initialize(num_genes,prior_g,spp);
    in.clear();
    in.close();
    cout << "still populating" << endl;
    in.open(file);
    genome_index=0;
    while(in.peek()!=EOF)
    {
        //cout << g[genome_index].acc_genome << "\t" << g[genome_index].spp_information << endl;
        for(i=0;i<g[genome_index].num_proteins;i++)
        {
            in.getline(text,100,'\t');  //strip off the accession number
            in.getline(text,100,'\t');  //get protein id number
            g[genome_index].add_protein(i,text);
            in.getline(text,100,'\n');
        }
        genome_index++;
        //cout << "Looping populating" << endl;
    }
    cout << "Text " << text << endl;
    cout << "Done populating" << endl;
    in.clear();
    in.close();
    cout << "Text " << text << endl;
    delete text;
    delete prior_g;
    return true;
}

bool super_file::read_BLAST(char * file)
{
    cout << "Reading BLAST" << endl;
    ifstream in;
    in.open(file);

    if(!in.is_open()) {cout << "Cannot find file.\n"; return false;}

    char node[1000];
    char *hit=new char [1000];
    char junk[50];
    float evalue;
    int genome_num, gene_num;

    char * prot_acc_num;
    char * prot_function;
    char * spp;

    while(in.peek()!=EOF)
    {
        in.getline(node,1000,'\t');
        in.getline(hit,1000,'\t');
        in>>evalue;
        in.getline(junk,50,'\n');
        prot_acc_num=NULL; prot_function=NULL; spp=NULL;
        if(parse_hit(hit,prot_acc_num,prot_function,spp))
        {
            //update count
            if(is_in(prot_acc_num, genome_num, gene_num))
                g[genome_num].increment_hits_counter(gene_num);

        }
        delete prot_acc_num; delete prot_function; delete spp;

    }
    in.clear();
    in.close();

    delete hit;
    delete prot_acc_num;
    delete prot_function;
    delete spp;

    return true;
}

bool super_file::read_BLAST(char * file, float evalue_threshold)
{
    cout << "Reading BLAST 2-param" << endl;

    ifstream in;
    in.open(file);
    cout << "seemed to have opened file " << file << endl;
    if(!in.is_open()) {cout << "Cannot find file.\n"; return false;}

    char node[1000];
    char *hit=new char [1000];
    char junk[50];
    float evalue;
    int genome_num, gene_num;

    char * prot_acc_num;
    char * prot_function;
    char * spp;
    cout << "Reading file " << endl;

    while(in.peek()!=EOF)
    {
        in.getline(node,1000,'\t');
        in.getline(hit,1000,'\t');
        in>>evalue;
        in.getline(junk,50,'\n');
        if(evalue<=evalue_threshold)
        {
            prot_acc_num=NULL; prot_function=NULL; spp=NULL;
            if(lambda)
            {
	      if(parse_hit(hit,prot_acc_num))
                {
                    //update count
                    if(is_in(prot_acc_num, genome_num, gene_num))
                        g[genome_num].increment_hits_counter(gene_num);
                }
                delete prot_acc_num;
            }
            if(pauda)
            {
                if(parse_hit(hit,prot_acc_num,prot_function,spp))
                {
                    //update count
                    if(is_in(prot_acc_num, genome_num, gene_num))
                        g[genome_num].increment_hits_counter(gene_num);
                }
                delete prot_acc_num; delete prot_function; delete spp;
            }
            if(other)
            {

            }
        }
    }
    cout << "Finished file" << endl;
    in.clear();
    in.close();
    delete hit;
 
    cout << "Finished reading BLAST" << endl;
    return true;
}


bool super_file::is_in(char * acc, int & genome_num, int & gene_num)
{
    int i,j;
    for(i=0;i<num_genomes;i++)
    {
        for(j=0;j<g[i].num_proteins;j++)
        {
            if(strcmp(g[i].proteins[j],acc)==0)
            {
                genome_num=i; gene_num=j;
                return true;
            }
        }
    }
    genome_num=gene_num=0;
    return false;
}

bool super_file::write_out_forHitViz(char * output_file, char * output_root)
{
    cout << "Writing out for HitViz" << endl;
    int i,j,n;
    ofstream out;
    char file_name[300];

    //output #hits by protein
    sprintf(file_name,"%sHV_%s", output_root, output_file);
    cout << "File out: " << file_name << endl;
    out.open(file_name);
    for(i=0;i<num_genomes;i++)
    {
        n=g[i].get_number_of_hits_to_genome();
        if(n>0)
        {
            out << g[i].acc_genome << "|" << g[i].spp_information << "|" << endl;
            for(j=0;j<g[i].num_proteins;j++)
            {
                out << g[i].proteins[j] << "|" << g[i].counts[j] << "|" << endl;
            }
            out << "*" << endl;
        }
    }
    out.clear();
    out.close();

    return true;
}

bool super_file::write_out_forKrona(char * output_file, char * output_root)
{
    cout << "Writing out for Krona" << endl;
    int i;
    ofstream out;
    char file_name[300];
    //output coverage by genome
    sprintf(file_name,"%sKronaStats_%s",output_root, output_file);
    cout << "File out: " << file_name << endl;
    out.open(file_name);
    out << "Accession Number\tSpp\tNumber of Proteins in Genome\tNumber of Hits to Genome\t% of Proteins Hit\n";
    out << "Number of genomes " << num_genomes << endl;
    for(i=0;i<num_genomes;i++)
    {
        out << g[i].acc_genome << "\t" << g[i].spp_information << "\t" << g[i].num_proteins << "\t";
        out << g[i].get_number_of_hits_to_genome() << "\t" << g[i].get_percentage_of_proteins_hit() << endl;
    }
    out.clear();
    out.close();
    cout << "Finished Krona" << endl;
    return true;
}

bool super_file::write_out_summary_statistics(char * output_file, char * output_root)
{
    cout << "Writing out summary statistics" << endl;
    int i,j,n;
    ofstream out;
    char file_name[300];

    //output coverage by genome
    sprintf(file_name,"%scoverage_%s",output_root, output_file);
    cout << "File out: " << file_name << endl;

    out.open(file_name);
    out << "Accession Number\tSpp\tNumber of Proteins in Genome\tNumber of Hits to Genome\t% of Proteins Hit\n";
    out << "Number of genomes " << num_genomes << endl;
    for(i=0;i<num_genomes;i++)
    {
        out << g[i].acc_genome << "\t" << g[i].spp_information << "\t" << g[i].num_proteins << "\t";
        out << g[i].get_number_of_hits_to_genome() << "\t" << g[i].get_percentage_of_proteins_hit() << endl;
    }
    out.clear();
    out.close();

    //output #hits by protein
    sprintf(file_name,"%shits_by_protein_%s",output_root, output_file);
    out.open(file_name);
    for(i=0;i<num_genomes;i++)
    {
        n=g[i].get_number_of_hits_to_genome();
        if(n>0)
        {
            out << g[i].acc_genome << "\t" << g[i].spp_information << "\t" << n << endl;
            out << "Protein Acc #\tNumber of Hits to Protein\n";
            for(j=0;j<g[i].num_proteins;j++)
            {
                if(g[i].counts[j]>0) out << g[i].proteins[j] << "\t" << g[i].counts[j] << endl;
            }
            out << endl;
        }
    }
    out.clear();
    out.close();
    cout << "Finished summary stats" << endl;
    return true;
}

bool super_file::reset_counts()
{
   int i;
   for(i=0;i<num_genomes;i++)
        g[i].initialize_counts();
   return true;
}

super_file::~super_file()
{
   for(int i=0;i<num_genomes;i++)
     g[i].~genome();
}


/** PARSE HIT DATA FOR LAMBDA BLASTS **/
bool super_file::parse_hit(char *& hit, char *& prot_acc_num)
{
    int i,j,start_p,end_p;
    char * pch;
    bool special_case=false;

    if(hit[1]!='i') return false;
    pch=strchr(hit,'|');
    if(pch==NULL) return false;
    pch=strchr(pch+1,'|');  //get second pipe

    //get accession
    pch=strchr(pch+1,'|');  //get third pipe
    start_p=pch-hit;
    pch=strchr(pch+1,'|');  //get fourth pipe
    end_p=pch-hit;
    j=end_p-start_p;
    prot_acc_num=new char [j+1];
    for(i=0;i<j-1;i++) prot_acc_num[i]=hit[start_p+1+i];
    prot_acc_num[i]='\0';
    return true;
}

/** PARSE HIT DATA FOR PAUDA BLASTS **/
bool super_file::parse_hit(char *& hit, char *& prot_acc_num, char *& prot_function, char *& spp)
{
    int i,j,start_p,end_p;
    char * pch;

    pch=strstr(hit,"gi|");
    if(pch==NULL)  return false;
    else
    {
        pch=strchr(hit,'|');
        if(pch==NULL) return false;
        pch=strchr(pch+1,'|');  //get second pipe
        if(pch==NULL) return false;
        //get accession
        pch=strchr(pch+1,'|');  //get third pipe
        if(pch==NULL) return false;
        start_p=pch-hit;
        pch=strchr(pch+1,'|');  //get fourth pipe
        if(pch==NULL) return false;
        end_p=pch-hit;
        j=end_p-start_p;
        prot_acc_num=new char [j+1];
        for(i=0;i<j-1;i++) prot_acc_num[i]=hit[start_p+1+i];
        prot_acc_num[i]='\0';

        //get function
        pch=strchr(pch+1,'[');  //get bracket
        if(pch==NULL) return false;

        start_p=end_p+1;
        end_p=pch-hit;
        j=end_p-start_p+1;
        prot_function=new char [j+1];
        while(hit[start_p]==' ') {start_p++; j--;}
        for(i=0;i<j-1;i++) prot_function[i]=hit[start_p+i];
        if(prot_function[i-1]==' ') i--;
        prot_function[i]='\0';

        //get species
        pch=strchr(pch+1,']');  //get bracket
        start_p=end_p+1;
        end_p=pch-hit;
        j=end_p-start_p+1;
        spp=new char[j+1];
        for(i=0;i<j-1;i++) spp[i]=hit[start_p+i];
        spp[i]='\0';

        return true;
    }
}



/***************************************************************************************************************
function creates a list of genome acc number, protein acc number, virus species information. separated by tabs.
***************************************************************************************************************/

int make_protein_list_file(char * file_of_files)
{
    cout << "Making protein list" << endl;
    cout.flush();
    int i,ii,j,k,start_p,end_p;
    char ** files;
    int num_faa_files=0;
    char g_acc[2000];
    char prot_acc[2000];
    char g_name[2000];
    char line[2000];
    ifstream in;
    in.open(file_of_files);
    while(in.peek()!=EOF)
    {
        in.getline(line,2000,'\n');
        num_faa_files++;
    }
    in.clear();
    in.close();

    files=new char * [num_faa_files];
    for(i=0;i<num_faa_files;i++) files[i]=new char [500];
    in.open(file_of_files);
    for(i=0;i<num_faa_files;i++)
    {
        in.getline(files[i],500,'\n');
    }
    in.clear();
    in.close();

    ofstream out;
    out.open("genome_proteins_list_full2.txt");
    char * pch;

    for(i=0;i<num_faa_files;i++)
    {
        in.open(files[i]);

        //parse out file name
        pch=strrchr(files[i],'/');
        if(pch==NULL) j=0;
        else j=pch-files[i];
        for(k=j; k<strlen(files[i]); k++)
        {
            g_acc[k-j]=files[i][k];
            if(g_acc[k-j]=='.') {g_acc[k-j]='\0'; k=strlen(files[i]);}
        }
        g_acc[k-j]='\0';

        while(in.peek()!=EOF)
        {
            in.getline(line,2000,'\n');
            if(line[0]=='>')
            {
                pch=strrchr(line,'|');      //first pipe
                if(pch==NULL)  //not in standard output
                    out << g_acc << "\t" << line << "\t" << line << endl;
                else
                {
                    //parse out protein acc number
                    pch=strchr(line,'|');   //first pipe
                    pch=strchr(pch+1,'|');  //second pipe
                    pch=strchr(pch+1,'|');  //third pipe
                    start_p=pch-line;
                    pch=strchr(pch+1,'|');  //get fourth pipe
                    end_p=pch-line;
                    j=end_p-start_p;
                    //prot_acc=new char [j+1];
                    for(ii=0;ii<j-1;ii++)
                        prot_acc[ii]= line[start_p+1+ii];
                    prot_acc[ii]='\0';

                    //parse out species information
                    pch=strrchr(line,'[');
                    if(pch==NULL)   //not in standard output
                        sprintf(g_name,"%s",line);
                    else
                    {
                        j=pch-line+1;
                        pch=strrchr(line,']');
                        if(pch==NULL)   //not in standard output
                            sprintf(g_name,"%s",line);
                        else
                        {
                            k=pch-line;
                            for(ii=j;ii<k;ii++) g_name[ii-j]=line[ii];
                            g_name[ii-j]='\0';
                        }
                    }
                    out << g_acc << "\t" << prot_acc << "\t" << g_name << endl;
                }
            }
        }
        in.clear();
        in.close();
    }
    out.clear();
    out.close();
    for(i=0;i<num_faa_files;i++) delete files[i];
    delete files;
    return num_faa_files;
}




/**************************************************************************************
argument 0= BLASTx (or BLASTx-like) tool
argument 1= PAUDA/LAMBDA indicator
argument 2= evalue threshold
argument 3= file which is a list of the path and file names of all of the FASTA format
            faa genome files used for comparison (db) 
argument 4= input parsed ORF file from Pauda or Lambda
argument 5= output file name root for files to be named
***************************************************************************************/
int main(int argc, char * argv[]) {

    cout << "Starting to run... hopefully" << endl;
    cout.flush();
    int i;

     //VARIABLES FOR TESTING WITHOUT USING THE COMMAND LINE
    
    //    argv[1]="LAMBDA";
    //argv[2]="0.0001";
    //argv[3]="/home/lsb456/virusLand2/faa/faaList.txt";
    //argv[4]="/home/lsb456/virusLand2/run2/lambda/v1S1Xv1S2_ORF_chopped.m8";
    //argv[5]="myTest";//"/home/lsb456/virusLand2/run2/stats/";
    

    //determine the means in which virus reads were classified
    char method=argv[1][0];
    if(method=='P') {cout << "It's PAUDA" << endl; pauda=true; lambda=false; other=false;}
    if(method=='L') {cout << "It's LAMBDA" << endl; pauda=false; lambda=true; other=false;}
    if(method=='O') {cout << "It's OTHER" << endl; pauda=false; lambda=false; other=true;}


//    if(argv[1]=="PAUDA") { cout << "It's PAUDA" << endl; pauda=true; lambda=false; other=false;}
//    if(argv[1]=="LAMBDA") {cout << "It's LAMBDA" << endl; pauda=false; lambda=true; other=false;}
//    if(argv[1]=="OTHER") {cout << "It's OTHER" << endl; pauda=false; lambda=false; other=true;}

    //char * pauda_output_file = "/home/lsb456/virusLand/paudaParsed/v1S1Xv1S4.parsd";
    float threshold = atof(argv[2]);
cout << threshold << endl << endl;
    cout << "*****" << endl;
    cout << "Reading FAA list file: " << argv[3] << endl;
    int num_faa_file= make_protein_list_file(argv[3]);
    cout << "Pauda file" << argv[4] << endl;

     if(num_faa_file>0)
      {
        cout << "Doing something! I should have executed 'make_protein_list_file' by now" << endl;
        cout.flush();
       	super_file s(num_faa_file,"genome_proteins_list_full2.txt");

        s.read_BLAST(argv[4],threshold);
	cout << "done with read_BLAST" << endl; cout.flush();

        //for HitViz
        s.write_out_forHitViz(argv[1],argv[5]);

        //for Krona
	s.write_out_forKrona(argv[1], argv[5]);

        //for Statistics Files
	s.write_out_summary_statistics(argv[1], argv[5]);
    }
    else cout << "Could Not Read Library FAA FILES.\n";

    return 1;
}
