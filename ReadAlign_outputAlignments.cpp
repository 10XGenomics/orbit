#include "ReadAlign.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"

void ReadAlign::outputAlignments() {
    
    outBAMbytes=0;
    outBAMbytes1=0;

    bool mateMapped[2]={false,false};
    
    if (P->outFilterBySJoutStage<=1) {//no chimeric output for stage=2
        chimericDetection();
    };
    
    if ( nW==0 ) {//no good windows
        statsRA.unmappedOther++;
        unmapType=0;
    } else if ( (trBest->maxScore < P->outFilterScoreMin) || (trBest->maxScore < (intScore) (P->outFilterScoreMinOverLread*(Lread-1))) \
              || (trBest->nMatch < P->outFilterMatchNmin)  || (trBest->nMatch < (uint) (P->outFilterMatchNminOverLread*(Lread-1))) ) {//too short
        statsRA.unmappedShort++;
        unmapType=1;
    } else if ( (trBest->nMM > outFilterMismatchNmaxTotal) || (double(trBest->nMM)/double(trBest->rLength)>P->outFilterMismatchNoverLmax) ) {//too many mismatches
        statsRA.unmappedMismatch++;
        unmapType=2;
    } else if (nTr > P->outFilterMultimapNmax){//too multi
        statsRA.unmappedMulti++;
        unmapType=3;
    } else {//output transcripts

        outFilterPassed=true;
        if (P->outFilterBySJoutStage==1) {//filtering by SJout
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
                statsRA.readN--;
                statsRA.readBases -= readLength[0]+readLength[1];
                
//                 if (P->runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexOutFilterBySJout);
                for (uint im=0;im<P->readNmates;im++) {
                   chunkOutFilterBySJoutFiles[im] << readNameMates[im] << "/" <<im+1 <<" "<< iReadAll <<" "<< readFilter <<" "<< readFilesIndex;
                   chunkOutFilterBySJoutFiles[im] <<"\n";
                   chunkOutFilterBySJoutFiles[im] << Read0[im] <<"\n";
                    if (readFileType==2) {//fastq
                        chunkOutFilterBySJoutFiles[im] << "+\n";
                        chunkOutFilterBySJoutFiles[im] << Qual0[im] <<"\n";
                    };
                };
//                 if (P->runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexOutFilterBySJout);  
            };
        };

        if (P->outSJfilterReads=="All" || nTr==1) {
            uint sjReadStartN=chunkOutSJ1->N;        
            for (uint iTr=0;iTr<nTr;iTr++) {//write all transcripts
                outputTranscriptSJ (*(trMult[iTr]), nTr, chunkOutSJ1, sjReadStartN);            
            };
        };        

        if (outFilterPassed) {
            if (nTr>1) {//multimappers
                statsRA.mappedReadsM++;
                unmapType=-1;
            } else if (nTr==1) {//unique mappers
                statsRA.mappedReadsU++;
                statsRA.transcriptStats(*(trMult[0]),Lread);
                unmapType=-2;
            } else {//cannot be
                ostringstream errOut;
                errOut  << "EXITING because of a BUG: nTr=0 in outputAlignments.cpp";
                exitWithError(errOut.str(), std::cerr, P->inOut->logMain, EXIT_CODE_BUG, *P);                    
            };            
            
            for (uint iTr=0;iTr<nTr;iTr++) {//write all transcripts
                if (P->outBAMunsorted || P->outBAMcoord) {//BAM output
                    alignBAM(*(trMult[iTr]), nTr, iTr, P->chrStart[trMult[iTr]->Chr], (uint) -1, (uint) -1, 0, -1, NULL, P->outSAMattrOrder);
                    for (uint imate=0; imate<P->readNmates; imate++) {//output each mate
                        if (P->outBAMunsorted) outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                        if (P->outBAMcoord)    outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], P->chrStart[trMult[iTr]->Chr], (iReadAll<<32) | (iTr<<8) | trMult[iTr]->exons[0][EX_iFrag] );                        
                    };
                };
            
                if (P->outSAMbool){//SAM output
                    outBAMbytes+=outputTranscriptSAM(*(trMult[iTr]), nTr, iTr, (uint) -1, (uint) -1, 0, -1, NULL, outSAMstream);
                };
            };
                        
            if (P->outSJfilterReads=="All" || nTr==1) {
                uint sjReadStartN=chunkOutSJ->N;        
                for (uint iTr=0;iTr<nTr;iTr++) {//write all transcripts
                    outputTranscriptSJ (*(trMult[iTr]), nTr, chunkOutSJ, sjReadStartN);            
                };
            };
            mateMapped[trBest->exons[0][EX_iFrag]]=true;
            mateMapped[trBest->exons[trBest->nExons-1][EX_iFrag]]=true;        
            if (P->readNmates>1 && !(mateMapped[0] && mateMapped[1]) ) unmapType=4;
            
            if ( P->quant.trSAM.yes ) {
                quantTranscriptome(chunkTr, nTr, trMult,  alignTrAll);
            };
            if ( P->quant.geCount.yes ) {
                chunkTr->geneCountsAddAlign(nTr, trMult);
            };           
        };
    };

    if (unmapType>=0 && P->outSAMunmapped=="Within") {//unmapped read, at least one mate
        if (P->outBAMunsorted || P->outBAMcoord || P->quant.trSAM.yes) {//BAM output
            alignBAM(*trBest, 0, 0, P->chrStart[trBest->Chr], (uint) -1, (uint) -1, 0,  unmapType, mateMapped, P->outSAMattrOrder);
            for (uint imate=0; imate<P->readNmates; imate++) {//output each mate
                if (P->outBAMunsorted) outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                //TODO clean for single-end alignments of PE reads
                if ( P->quant.trSAM.yes && unmapType!=4) outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                if (P->outBAMcoord)    outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], 0, iReadAll);                                        
            };
        };
        if (P->outSAMbool) {
            outBAMbytes+= outputTranscriptSAM(*trBest, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);        
        };        
    };
    if (unmapType>=0 && P->outReadsUnmapped=="Fastx" ){//output to fasta/q files
           for (uint im=0;im<P->readNmates;im++) {
               chunkOutUnmappedReadsStream[im] << readNameMates[im] << "/" <<im+1;
               if (P->readNmates>1) chunkOutUnmappedReadsStream[im] <<"\t"<< int(mateMapped[0]) <<  int(mateMapped[1]);
               chunkOutUnmappedReadsStream[im] <<"\n";
               chunkOutUnmappedReadsStream[im] << Read0[im] <<"\n";
                if (readFileType==2) {//fastq
                    chunkOutUnmappedReadsStream[im] << "+\n";
                    chunkOutUnmappedReadsStream[im] << Qual0[im] <<"\n";
                };
           };
    }; 
};



