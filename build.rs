extern crate bindgen;
extern crate cc;
extern crate skeptic;
use std::process::Command;
use std::path::Path;

fn main() {

    let status = Command::new("make").current_dir("STAR/source")
        .arg("liborbit.a")
        .status().expect("Failed to build STAR");
    assert!(status.success());
    

    println!("cargo:rustc-link-search=rust-utils-10x/orbit/STAR/source");
    println!("cargo:rustc-link-lib=static=orbit");
    //println!("cargo:rustc-link-search=/mnt/build/toolchains/versions/2.2.6/lib");
    /*
    let path = Path::new("STAR/source");
    let objects = vec!["ParametersSolo.o", "SoloRead.o", "SoloRead_record.o", "SoloReadBarcode.o", "SoloReadBarcode_getCBandUMI.o", "SoloReadFeature.o", "SoloReadFeature_record.o", "SoloReadFeature_inputRecords.o", "Solo.o", "SoloFeature.o", "SoloFeature_collapseUMI.o", "SoloFeature_outputResults.o", "SoloFeature_processRecords.o", "ReadAlign_outputAlignments.o", "ReadAlign.o", "SharedMemory.o", "PackedArray.o", "SuffixArrayFuns.o", "Parameters.o", "InOutStreams.o", "SequenceFuns.o", "Genome.o", "Stats.o", "Transcript.o", "Transcript_alignScore.o", "Transcript_generateCigarP.o", "Chain.o", "Transcript_variationAdjust.o", "Variation.o", "ReadAlign_waspMap.o", "ReadAlign_storeAligns.o", "ReadAlign_stitchPieces.o", "ReadAlign_multMapSelect.o", "ReadAlign_mapOneRead.o", "readLoad.o", "ReadAlignChunk.o", "OutSJ.o", "outputSJ.o", "blocksOverlap.o", "ThreadControl.o", "sysRemoveDir.o", "ReadAlign_maxMappableLength2strands.o", "binarySearch2.o", "ReadAlign_outputTranscriptSAM.o", "ReadAlign_outputTranscriptSJ.o", "ReadAlign_outputTranscriptCIGARp.o", "ReadAlign_createExtendWindowsWithAlign.o", "ReadAlign_assignAlignToWindow.o", "ReadAlign_oneRead.o", "ReadAlign_stitchWindowSeeds.o", "ReadAlign_peOverlapMergeMap.o", "ReadAlign_mappedFilter.o", "ReadAlign_chimericDetection.o", "ReadAlign_chimericDetectionOld.o", "ReadAlign_chimericDetectionOldOutput.o", "ChimericDetection.o", "ChimericDetection_chimericDetectionMult.o", "ReadAlign_chimericDetectionPEmerged.o", "ChimericAlign.o", "ChimericSegment.o", "ReadAlign_calcCIGAR.o", "ChimericAlign_chimericJunctionOutput.o", "ChimericAlign_chimericStitching.o", "stitchWindowAligns.o", "extendAlign.o", "stitchAlignToTranscript.o", "alignSmithWaterman.o", "Genome_genomeGenerate.o", "genomeParametersWrite.o", "genomeScanFastaFiles.o", "genomeSAindex.o", "Genome_insertSequences.o", "insertSeqSA.o", "funCompareUintAndSuffixes.o", "funCompareUintAndSuffixesMemcmp.o", "TimeFunctions.o", "ErrorWarning.o", "loadGTF.o", "streamFuns.o", "stringSubstituteAll.o", "Transcriptome.o", "Transcriptome_quantAlign.o", "Transcriptome_geneFullAlignOverlap.o", "ReadAlign_quantTranscriptome.o", "Quantifications.o", "Transcriptome_geneCountsAddAlign.o", "sjdbLoadFromFiles.o", "sjdbLoadFromStream.o", "sjdbPrepare.o", "sjdbBuildIndex.o", "sjdbInsertJunctions.o", "Parameters_openReadsFiles.o", "Parameters_closeReadsFiles.o", "Parameters_readSAMheader.o", "serviceFuns.o", "GlobalVariables.o", "BAMoutput.o", "BAMfunctions.o", "ReadAlign_alignBAM.o", "BAMbinSortByCoordinate.o", "BAMbinSortUnmapped.o", "orbit.o"];
    let cpps = objects.into_iter()
        .map(|o| path.join(o).with_extension("cpp"));

    cc::Build::new()
        .cpp(true)
        .warnings(false)
        .pic(true)
        .include("STAR/source")
        .files(cpps)
        .flag("-O3")
        .flag("-g")
        .flag("-pipe")
        .flag("-std=c++11")
        .flag("-fPIC")
        .flag("-DCOMPILATION_TIME_PLACE=\"/dev/null\"")
        .compile("orbit");
    */
    println!("cargo:rustc-link-lib=dylib=c++");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("wrapper.h")
        .whitelist_function("align_read")
        .whitelist_function("align_read_pair")
        .whitelist_function("init_aligner_clone")
        .whitelist_function("init_aligner")
        .whitelist_function("destroy_aligner")
        .whitelist_type("Aligner")
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file("src/bindings.rs")
        .expect("Couldn't write bindings!");

}
