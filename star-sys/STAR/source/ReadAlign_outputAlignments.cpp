#include "ReadAlign.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"

const char* ReadAlign::outputAlignments() {
    outBAMbytes=0;

    bool mateMapped[2]={false,false};

    set<uint32> readGeneFull={},readGene={};
    vector<uint32> readTranscripts={};
    vector<int32> readGeneExon={};

    std::stringstream stream;

    outFilterPassed=true;//only false if the alignment is held for outFilterBySJoutStage
    if (unmapType==-1) {//output transcripts
        if (P.outFilterBySJoutStage==1) {//filtering by SJout
            for (uint iTr=0;iTr<nTr;iTr++) {//check transcript for unannotated junctions
                for (uint iex=0;iex<trMult[iTr]->nExons-1;iex++) {//check all junctions
                    if (trMult[iTr]->canonSJ[iex]>=0 && trMult[iTr]->sjAnnot[iex]==0) {
                        outFilterPassed=false;
                        break;
                    };
                };
                if (!outFilterPassed) break;
            };
            if (!outFilterPassed) {//this read is held for further filtering BySJout, record fastq
                unmapType=-3; //the read is not conisddred unmapped
                for (uint im=0;im<readNmates;im++) {
                   chunkOutFilterBySJoutFiles[im] << readNameMates[im] <<" "<< iReadAll <<" "<< readFilter <<" "<< readFilesIndex;
                   if (!readNameExtra[im].empty())
                       chunkOutFilterBySJoutFiles[im]<<" "<< readNameExtra[im];
                   chunkOutFilterBySJoutFiles[im] <<"\n";
                   chunkOutFilterBySJoutFiles[im] << Read0[im] <<"\n";
                    if (readFileType==2) {//fastq
                        chunkOutFilterBySJoutFiles[im] << "+\n";
                        chunkOutFilterBySJoutFiles[im] << Qual0[im] <<"\n";
                    };
                };
            };
        };

        /*
        if (P.outSJfilterReads=="All" || nTr==1) {
            OutSJ *chunkOutSJ1=new OutSJ (P.limitOutSJcollapsed, P, mapGen);
            uint sjReadStartN=chunkOutSJ1->N;
            for (uint iTr=0;iTr<nTr;iTr++) {//report SJs for all transcripts
                //printf("outsj stuff\n");
                outputTranscriptSJ (*(trMult[iTr]), nTr, chunkOutSJ1, sjReadStartN);
            };
            delete chunkOutSJ1;
        };
        */

        if (outFilterPassed) {
            uint nTrOut=nTr; //number of aligns to output
            bool outSAMfilterYes=true;
            if (P.outSAMfilter.yes) {
                if (P.outSAMfilter.KeepOnlyAddedReferences) {
                    for (uint itr=0;itr<nTr;itr++) {//check if transcripts map to chr other than added references
                        if (trMult[itr]->Chr<mapGen.genomeInsertChrIndFirst) {
                            outSAMfilterYes=false;
                            break;
                        };
                    };
                } else if (P.outSAMfilter.KeepAllAddedReferences) {
                    nTrOut=0;
                    for (uint itr=0;itr<nTr;itr++) {//check if transcripts map to chr other than added references
                        if (trMult[itr]->Chr>=mapGen.genomeInsertChrIndFirst) {
                            trMult[nTrOut]=trMult[itr];
                            trMult[nTrOut]->primaryFlag=false;
                            ++nTrOut;
                        };
                    };
                    if (nTrOut==0) {
                        outSAMfilterYes=false;
                    } else {
                        trMult[0]->primaryFlag=true;
                    };
                };
            };
            if (nTr>1) {//multimappers
                unmapType=-1;
            } else if (nTr==1) {//unique mappers
                unmapType=-2;
            } else {//cannot be
                ostringstream errOut;
                errOut  << "EXITING because of a BUG: nTr=0 in outputAlignments.cpp";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
            };

            nTrOut=min(P.outSAMmultNmax,nTrOut); //number of to write to SAM/BAM files
            //write to SAM/BAM
            //printf("nTrOut %llu\n", nTrOut);
            for (uint iTr=0;iTr<nTrOut;iTr++) {//write all transcripts
                //mateMapped1 = true if a mate is present in this transcript
                bool mateMapped1[2]={false,false};
                mateMapped1[trMult[iTr]->exons[0][EX_iFrag]]=true;
                mateMapped1[trMult[iTr]->exons[trMult[iTr]->nExons-1][EX_iFrag]]=true;

                if (P.outSAMbool && outSAMfilterYes) {//SAM output
                    //printf("samout\n");
                    outBAMbytes+=outputTranscriptSAM(*(trMult[iTr]), nTr, iTr, (uint) -1, (uint) -1, 0, -1, NULL, &stream);
                    if (P.outSAMunmapped.keepPairs && readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) ) {//keep pairs && paired reads && one of the mates not mapped in this transcript
                        //printf("samout no null\n");
                        outBAMbytes+= outputTranscriptSAM(*(trMult[iTr]), 0, 0, (uint) -1, (uint) -1, 0, 4, mateMapped1, &stream);
                    };
                };
            };

            mateMapped[trBest->exons[0][EX_iFrag]]=true;
            mateMapped[trBest->exons[trBest->nExons-1][EX_iFrag]]=true;

            if (readNmates>1 && !(mateMapped[0] && mateMapped[1]) ) {
                unmapType=4;
            };


            if (unmapType==4 && P.outSAMunmapped.yes) {//output unmapped end for single-end alignments
                if (P.outSAMbool && !P.outSAMunmapped.keepPairs && outSAMfilterYes) {
                    outBAMbytes+= outputTranscriptSAM(*trBest, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, &stream);
                };
            };

            /*
            if (P.outSJfilterReads=="All" || nTr==1) {
                chunkOutSJ=new OutSJ (P.limitOutSJcollapsed, P, mapGen);
                uint sjReadStartN=chunkOutSJ->N;
                for (uint iTr=0;iTr<nTr;iTr++) {//write all transcripts junctions
                    outputTranscriptSJ (*(trMult[iTr]), nTr, chunkOutSJ, sjReadStartN);
                };
                delete chunkOutSJ;
            };
            */

            //transcripts
            if ( P.quant.trSAM.yes ) {//NOTE: the transcripts are changed by this function (soft-clipping extended), cannot be reused
                quantTranscriptome(chunkTr, nTrOut, trMult,  alignTrAll.get(), readTranscripts, readGene);
            };

        };
    };



    if ( P.outSAMunmapped.within && unmapType>=0 && unmapType<4 ) {//output unmapped within && unmapped read && both mates unmapped
        if (P.outSAMbool) {//output SAM
            outBAMbytes+= outputTranscriptSAM(*trBest, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, &stream);
            //printf("how about here?\n");
        };
    };

    std::stringbuf * pbuf = stream.rdbuf();
    std::streamsize size = pbuf->pubseekoff(0,stream.end);
    pbuf->pubseekoff(0,stream.beg);  
    char *contents = (char*)malloc(size+1);
    pbuf->sgetn (contents,size);
    contents[size] = '\0';

    return contents;
};



