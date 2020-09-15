// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

use std::env;
use std::path::Path;

fn libcxx() -> &'static str {
    match env::var("CXX") {
        Ok(cxx) => match Path::new(&cxx).file_name().unwrap().to_str().unwrap() {
            s if s.contains("clang++") => "c++",
            s if s.contains("g++") => "stdc++",
            s => panic!("unknown compiler: {}", s),
        },
        Err(_) => match env::var("TARGET") {
            Ok(ref s) if s.contains("darwin") => "c++",
            Ok(ref s) if s.contains("linux") => "stdc++",
            Ok(ref s) => panic!("unknown target: {}", s),
            Err(_) => panic!("TARGET is undefined"),
        },
    }
}

const FILES: &[&str] = &[
    "STAR/source/orbit.cpp",
    "STAR/source/InOutStreams.cpp",
    "STAR/source/Parameters.cpp",
    "STAR/source/ParametersSolo.cpp",
    "STAR/source/ParametersChimeric_initialize.cpp",
    "STAR/source/PackedArray.cpp",
    "STAR/source/MmapArray.cpp",
    "STAR/source/SequenceFuns.cpp",
    "STAR/source/Genome.cpp",
    "STAR/source/Genome_insertSequences.cpp",
    "STAR/source/Genome_genomeGenerate.cpp",
    "STAR/source/streamFuns.cpp",
    "STAR/source/genomeScanFastaFiles.cpp",
    "STAR/source/TimeFunctions.cpp",
    "STAR/source/insertSeqSA.cpp",
    "STAR/source/ReadAlign.cpp",
    "STAR/source/Transcript.cpp",
    "STAR/source/Transcriptome_quantAlign.cpp",
    "STAR/source/funCompareUintAndSuffixesMemcmp.cpp",
    "STAR/source/genomeSAindex.cpp",
    "STAR/source/ReadAlign_outputAlignments.cpp",
    "STAR/source/ReadAlign_outputTranscriptSAM.cpp",
    "STAR/source/ReadAlign_quantTranscriptome.cpp",
    "STAR/source/ReadAlign_calcCIGAR.cpp",
    "STAR/source/ReadAlign_storeAligns.cpp",
    "STAR/source/SuffixArrayFuns.cpp",
    "STAR/source/ReadAlign_oneRead.cpp",
    "STAR/source/ReadAlign_mapOneRead.cpp",
    "STAR/source/ReadAlign_stitchPieces.cpp",
    "STAR/source/ReadAlign_mappedFilter.cpp",
    "STAR/source/ReadAlign_maxMappableLength2strands.cpp",
    "STAR/source/ReadAlign_assignAlignToWindow.cpp",
    "STAR/source/ReadAlign_createExtendWindowsWithAlign.cpp",
    "STAR/source/ReadAlign_multMapSelect.cpp",
    "STAR/source/readLoad.cpp",
    "STAR/source/stitchWindowAligns.cpp",
    "STAR/source/extendAlign.cpp",
    "STAR/source/binarySearch2.cpp",
    "STAR/source/blocksOverlap.cpp",
    "STAR/source/stitchAlignToTranscript.cpp",
    "STAR/source/ErrorWarning.cpp",
];

// I used the following monstrosity to generate this list:
// head -62 ../build.rs | tail -40 | sed -e 's|\s\+"STAR\/||' | sed -e 's|",||' | xargs -I{} grep "#include" {} | sort | uniq | grep -v "<*>" | grep -v "^\/\/" | sed -e 's|^\s\+||' | sed -e 's|SAMTOOLS_BGZF_H|"htslib/htslib/bgzf.h"|' | sed -e 's|#include\ "|"STAR\/source\/|' | sed -e 's|"$|",|' | sort >> ../build.rs
const HEADERS: &[&str] = &[
    "STAR/source/alignSmithWaterman.h",
    "STAR/source/binarySearch2.h",
    "STAR/source/blocksOverlap.h",
    "STAR/source/ErrorWarning.h",
    "STAR/source/extendAlign.h",
    "STAR/source/funCompareUintAndSuffixes.h",
    "STAR/source/funCompareUintAndSuffixesMemcmp.h",
    "STAR/source/Genome.h",
    "STAR/source/genomeParametersWrite.h",
    "STAR/source/genomeSAindex.h",
    "STAR/source/genomeScanFastaFiles.h",
    "STAR/source/GlobalVariables.h",
    "STAR/source/htslib/htslib/bgzf.h",
    "STAR/source/IncludeDefine.h",
    "STAR/source/InOutStreams.h",
    "STAR/source/insertSeqSA.h",
    "STAR/source/loadGTF.h",
    "STAR/source/orbit.h",
    "STAR/source/OutSJ.h",
    "STAR/source/MmapArray.h",
    "STAR/source/PackedArray.h",
    "STAR/source/ParametersChimeric.h",
    "STAR/source/parametersDefault.xxd",
    "STAR/source/Parameters.h",
    "STAR/source/ParametersSolo.h",
    "STAR/source/ReadAlign.h",
    "STAR/source/readLoad.h",
    "STAR/source/SequenceFuns.h",
    "STAR/source/serviceFuns.cpp",
    "STAR/source/sjAlignSplit.cpp",
    "STAR/source/SjdbClass.h",
    "STAR/source/sjdbInsertJunctions.h",
    "STAR/source/sjdbLoadFromFiles.h",
    "STAR/source/sjdbPrepare.h",
    "STAR/source/sortSuffixesBucket.h",
    "STAR/source/Stats.h",
    "STAR/source/stitchWindowAligns.h",
    "STAR/source/streamFuns.h",
    "STAR/source/stringSubstituteAll.h",
    "STAR/source/SuffixArrayFuns.h",
    "STAR/source/sysRemoveDir.h",
    "STAR/source/TimeFunctions.h",
    "STAR/source/Transcript.h",
    "STAR/source/Transcriptome.h",
];

fn main() {
    // do not rebuild constitutively
    for file in FILES {
        println!("cargo:rerun-if-changed={}", file);
    }
    for file in HEADERS {
        println!("cargo:rerun-if-changed={}", file);
    }
    cc::Build::new()
        .cpp(true)
        .cpp_link_stdlib(Some(libcxx()))
        .define("COMPILATION_TIME_PLACE", "\"build.rs\"")
        .files(FILES)
        .flag("-std=c++11")
        .flag("-Wall")
        .flag("-Wextra")
        .flag("-Werror")
        .compile("orbit");
}
