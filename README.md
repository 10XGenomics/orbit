# Orbit

Orbit is a wrapper around STAR written in Rust.  Given a built STAR genome
index, Orbit can be used to align individual reads or pairs of reads and get
BAM-format alignment records as output.


The core alignment library is contained in src/lib.rs, along with a number of
tests which demonstrate its usage.  Build settings are contained in build.rs,
and bindgen is used to link the STAR API written in C and allow access to its
functions in src/bindings.rs.


## Notes on bindgen.

We took bindgen out of the build process to reduce build time 
because the bindings are so small. If you need to update the bingings
you can use the commandline version of bindgen and check in the resulting
bindings.rs.

For posterity, here's the old `build.rs` snippet:

```rust
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
```