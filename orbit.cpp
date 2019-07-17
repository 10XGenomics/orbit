#include <sstream>
#include <string>

#include "Genome.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "Transcriptome.h"
#include "Variation.h"

class Aligner
{
    public:
        Parameters *p;
        ReadAlign *ra;
};

string align_read(Aligner a, char **Read1, char **Qual1, unsigned long long read_length)
{
    a.p->readNmates = 1;
    a.ra->readNmates = 1;
    a.ra->Read0 = Read1;
    a.ra->Qual0 = Qual1;
    a.ra->Lread = read_length;
    a.ra->readLength[0] = read_length;
    a.ra->readLength[1] = read_length;
    int readStatus = a.ra->oneRead();
    if(readStatus != 0)
    {
        return "";
    }
    string str = a.ra->outputAlignments();
    return str;
}

Aligner init(int argInN, char* argIn[])
{
    Parameters *P = (Parameters*)malloc(sizeof(Parameters));
    new(P) Parameters();
    P->inputParameters(argInN, argIn);

    Genome *mainGenome = (Genome*)malloc(sizeof(Genome));
    new(mainGenome) Genome(*P);
    mainGenome->genomeLoad();

    Transcriptome *mainTranscriptome = NULL;

    mainGenome->Var = new Variation(*P, mainGenome->chrStart, mainGenome->chrNameIndex);

    ReadAlign *RA = (ReadAlign*)malloc(sizeof(ReadAlign));
    new(RA) ReadAlign(*P, *mainGenome, mainTranscriptome, 0);

    Aligner res;
    res.ra = RA;
    res.p = P;
    return res;
}

int main()
{
    char* arr[] = {
            "STAR", "--genomeDir", "/mnt/opt/refdata_cellranger/GRCh38-3.0.0/star",
            "--outSAMmultNmax", "50",
            "--runThreadN", "1",
            "--readNameSeparator", "space",
            "--outSAMunmapped", "Within", "KeepPairs",
            "--outSAMtype", "SAM",
            "--outStd", "SAM",
            "--outSAMorder", "PairedKeepInputOrder",
    };
    int len = sizeof(arr) / sizeof(arr[0]);
    Aligner a = init(len, arr);

    std::string line;
    std::ifstream infile("1.fastq");
    int lineNum = 0;
    char** curRead = (char**)malloc(sizeof(char*));
    curRead[0] = (char*)malloc(500*sizeof(char));
    while(std::getline(infile, line))
    {
        //printf("line %s\n", line.c_str());
        //printf("lineNum %d\n", lineNum);
        if(lineNum%4 == 1)
        {
            strcpy(curRead[0], line.c_str());
        }
        else if(lineNum%4 == 3)
        {
            char** curQual = (char**)malloc(sizeof(char*));
            curQual[0] = (char*)malloc(500*sizeof(char));
            strcpy(curQual[0], line.c_str());
            printf("read = %s\n", curRead[0]);
            printf("qual = %s\n", curQual[0]);
            string bam_line = align_read(a, curRead, curQual, line.length());
            printf("%s\n", bam_line.c_str());
            free(curQual[0]);
            free(curQual);
        }
        lineNum++;
        if(lineNum == 100) break;
    }
    free(curRead[0]);
    free(curRead);
    return 0;
}
