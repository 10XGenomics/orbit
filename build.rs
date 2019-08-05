extern crate bindgen;
extern crate cc;
extern crate skeptic;
use std::process::Command;

fn main() {

    let status = Command::new("make").current_dir("STAR/source")
        .arg("liborbit.a")
        .status().expect("Failed to build STAR");
    assert!(status.success());
    
    println!("cargo:rustc-link-search=rust-utils-10x/orbit/STAR/source");
    println!("cargo:rustc-link-lib=static=orbit");
    //println!("cargo:rustc-link-search=/mnt/build/toolchains/versions/2.2.6/lib");
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
