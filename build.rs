extern crate bindgen;
extern crate cc;
extern crate skeptic;
extern crate fs_utils;

use fs_utils::copy::copy_directory;
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {

    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    if !out.join("STAR").exists() {
        copy_directory("STAR", &out).unwrap();
    }

    let mut cfg = cc::Build::new();
    cfg.warnings(false).static_flag(true).pic(true);

    let tool = cfg.get_compiler();
    let (cc_path, cflags_env) = (tool.path(), tool.cflags_env());
    println!("cflag: {:?}", cflags_env);
    let cc_cflags = cflags_env.to_string_lossy().replace("-O0", "");
    if Command::new("make")
        .current_dir(out.join("STAR/source"))
        .arg(format!("CC={}", cc_path.display()))
        .arg(format!("CFLAGS={}", cc_cflags))
        .arg("liborbit.a")
        .status()
        .unwrap()
        .success() != true
    {
        panic!("failed to build STAR");
    }
    
    let out_src = out.join("STAR").join("source");
    println!("cargo:rustc-link-search=native={}", out_src.display());
    println!("cargo:rustc-link-lib=static=orbit");
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
        .whitelist_function("init_star_ref")
        .whitelist_function("init_aligner_from_ref")
        .whitelist_function("destroy_aligner")
        .whitelist_function("destroy_ref")
        .whitelist_type("Aligner")
        .whitelist_type("StarRef")
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file("src/bindings.rs")
        .expect("Couldn't write bindings!");
}
