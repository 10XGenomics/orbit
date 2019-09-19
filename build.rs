// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

use fs_utils::copy::copy_directory;
use std::env;
use std::path::Path;
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
        .success()
        != true
    {
        panic!("failed to build STAR");
    }

    let out_src = out.join("STAR").join("source");
    println!("cargo:rustc-link-search=native={}", out_src.display());
    println!("cargo:rustc-link-lib=static=orbit");

    match env::var("CXX") {
        Ok(cxx) => {
            let cxx = Path::new(&cxx).file_name().unwrap().to_str().unwrap();
            if cxx.contains("clang++") {
                println!("cargo:rustc-link-lib=dylib=c++");
            } else if cxx.contains("g++") {
                println!("cargo:rustc-link-lib=dylib=stdc++");
            } else {
                panic!("unknown compiler: {}", cxx);
            }
        }
        Err(_) => {
            let target = env::var("TARGET").unwrap();
            if target.contains("darwin") {
                println!("cargo:rustc-link-lib=dylib=c++");
            } else if target.contains("linux") {
                println!("cargo:rustc-link-lib=dylib=stdc++");
            } else {
                panic!("unknown target: {}", target);
            }
        }
    }
}
