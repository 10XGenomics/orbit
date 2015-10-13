#include "Transcript.h"
#include "serviceFuns.h"

int Transcript::variationAdjust(Variation &Var, char *R)
{
    if (!Var.yes)
    {//no variation
        return 0;
    };
    
    int dScore=0;//chnage in the score
    
    //for each block, check whether it overlaps one or more SNPs
    for (uint ie=0; ie<nExons; ie++)
    {
        //binary search to find nearest SNP
        int32 isnp=binarySearch1b <uint> (exons[ie][EX_G], Var.snp.loci, Var.snp.N);
        if (isnp>0)
        {
            while (exons[ie][EX_G]+exons[ie][EX_L]>Var.snp.loci[isnp])
            {//these SNPs overlap the block
                snpInd.push_back(isnp); //record snp index
                snpLoci.push_back(Var.snp.loci[isnp]-Var.P->chrStart[Chr]);
                
                char ntR=R[exons[ie][EX_R]+Var.snp.loci[isnp]-exons[ie][EX_G]];//nt of the read in the SNP position, already trnasformed to + genome strand
                uint8 igt;
                for (igt=1; igt<3; igt++)
                {//1st or 2nd allele
                    if (Var.snp.nt[isnp][igt]==ntR)
                    {
                        break;
                    };
                };
                if (ntR == Var.snp.nt[isnp][0])
                {//mark snp that agrees with the reference
                    igt*=10;
                };
                snpGt.push_back(igt);
                
                if (igt<3 && ntR != Var.snp.nt[isnp][0])
                {//non-reference allele, correct nMM and score
                    dScore+=2;
                    nMM--;
                    nMatch++;
                };

                isnp++;
            };
        };
    };
    return dScore;
};
