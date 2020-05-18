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

fn main() {
    // do not rebuild constitutively
    for file in FILES {
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
