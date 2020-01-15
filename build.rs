// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

use fs_utils::copy::copy_directory;
use glob::glob;
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
    let out = env::var("OUT_DIR").unwrap();
    let out_source = format!("{}/source", out);
    match copy_directory("STAR/source", &out) {
        Ok(_) => (),
        Err(fs_utils::copy::Error::DestinationDirectoryExists(_)) => (),
        e => panic!(e),
    };

    let mut cfg = cc::Build::new();
    cfg.warnings(false).static_flag(true).pic(true);

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
}
