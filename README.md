# Orbit

Orbit is a wrapper around STAR written in Rust.  Given a built STAR genome
index, Orbit can be used to align individual reads or pairs of reads and get
BAM-format alignment records as output.


The core alignment library is contained in src/lib.rs, along with a number of
tests which demonstrate its usage.  Build settings are contained in build.rs,
and bindgen is used to link the STAR API written in C and allow access to its
functions in src/bindings.rs.
