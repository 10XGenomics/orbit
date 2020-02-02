// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

//use fs_utils::copy::copy_directory;
//use glob::glob;
use std::env;
use std::path::Path;
use std::process::Command;

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

fn main() {
    /*
    let out = env::var("OUT_DIR").unwrap();
    let out_source = format!("{}/source", out);

    match copy_directory("STAR/source", &out) {
        Ok(_) => (),
        Err(fs_utils::copy::Error::DestinationDirectoryExists(_)) => (),
        e => panic!(e),
    };
    */
    cc::Build::new()
        .cpp(true)
        .static_flag(true)
        .pic(true)
        .cpp_set_stdlib(Some(libcxx()))
        .define("COMPILATION_TIME_PLACE", "\"build.rs\"")
        .file("STAR/source/orbit.cpp")
        .file("STAR/source/InOutStreams.cpp")
        .file("STAR/source/Parameters.cpp")
        .file("STAR/source/ParametersSolo.cpp")
        .file("STAR/source/ParametersChimeric_initialize.cpp")
        .file("STAR/source/PackedArray.cpp")
        .file("STAR/source/SequenceFuns.cpp")
        .file("STAR/source/Genome.cpp")
        .file("STAR/source/Genome_insertSequences.cpp")
        .file("STAR/source/Genome_genomeGenerate.cpp")
        .file("STAR/source/streamFuns.cpp")
        .file("STAR/source/genomeScanFastaFiles.cpp")
        .file("STAR/source/TimeFunctions.cpp")
        .file("STAR/source/insertSeqSA.cpp")
        .file("STAR/source/ReadAlign.cpp")
        .file("STAR/source/Transcript.cpp")
        .file("STAR/source/Transcriptome_quantAlign.cpp")
        .file("STAR/source/funCompareUintAndSuffixesMemcmp.cpp")
        .file("STAR/source/genomeSAindex.cpp")
        .file("STAR/source/ReadAlign_outputAlignments.cpp")
        .file("STAR/source/ReadAlign_outputTranscriptSAM.cpp")
        .file("STAR/source/ReadAlign_quantTranscriptome.cpp")
        .file("STAR/source/ReadAlign_calcCIGAR.cpp")
        .file("STAR/source/ReadAlign_storeAligns.cpp")
        .file("STAR/source/SuffixArrayFuns.cpp")
        .file("STAR/source/ReadAlign_oneRead.cpp")
        .file("STAR/source/ReadAlign_mapOneRead.cpp")
        .file("STAR/source/ReadAlign_stitchPieces.cpp")
        .file("STAR/source/ReadAlign_mappedFilter.cpp")
        .file("STAR/source/ReadAlign_maxMappableLength2strands.cpp")
        .file("STAR/source/ReadAlign_assignAlignToWindow.cpp")
        .file("STAR/source/ReadAlign_createExtendWindowsWithAlign.cpp")
        .file("STAR/source/ReadAlign_multMapSelect.cpp")
        .file("STAR/source/readLoad.cpp")
        .file("STAR/source/stitchWindowAligns.cpp")
        .file("STAR/source/extendAlign.cpp")
        .file("STAR/source/binarySearch2.cpp")
        .file("STAR/source/blocksOverlap.cpp")
        .file("STAR/source/stitchAlignToTranscript.cpp")
        .file("STAR/source/ErrorWarning.cpp")
        .flag("-std=c++11")
        .flag("-Wall")
        .flag("-Wextra")
        .flag("-fPIC")
        .compile("orbit");

    /*
    let tool = cfg.get_compiler();
    let (cc_path, cflags_env) = (tool.path(), tool.cflags_env());
    println!("cflag: {:?}", cflags_env);
    let cc_cflags = cflags_env.to_string_lossy().replace("-O0", "");
    if Command::new("make")
        .current_dir(&out_source)
        .arg(format!("CC={}", cc_path.display()))
        .arg(format!("CFLAGS={}", &cc_cflags))
        .arg("liborbit.a")
        .status()
        .unwrap()
        .success()
        != true
    {
        panic!("failed to build STAR");
    }


    let from = format!("{}/liborbit.a", &out_source);
    let to = format!("{}/liborbit.a", &out);
    std::fs::rename(&from, &to).expect(&format!("{}: {}", from, to));

    std::fs::remove_dir_all(&out_source).expect(&out_source);

    println!("cargo:rustc-link-search=native={}", &out);
    println!("cargo:rustc-link-lib=static=orbit");
    println!("cargo:rustc-link-lib=dylib={}", libcxx());

    // so we don't rebuild every time
    for star_source in glob("STAR/source/**/*").expect("error in glob(STAR/source)") {
        let star_source = star_source
            .as_ref()
            .expect("error in glob(STAR/source)")
            .to_str()
            .expect("STAR source could not be made into valid unicode");
        println!("cargo:rerun-if-changed={}", star_source);
    }
    */
}
