//! Map Rust target features to the C++ compiler flags that enable the same instructions.

/// Rust target feature names and the flag that enables the same instructions in both
/// clang and gcc.
const FEATURE_FLAGS: [(&str, &str); 14] = [
    ("sse3", "-msse3"),
    ("ssse3", "-mssse3"),
    ("sse4.1", "-msse4.1"),
    ("sse4.2", "-msse4.2"),
    ("popcnt", "-mpopcnt"),
    ("cmpxchg16b", "-mcx16"),
    ("avx", "-mavx"),
    ("avx2", "-mavx2"),
    ("fma", "-mfma"),
    ("bmi1", "-mbmi"),
    ("bmi2", "-mbmi2"),
    ("lzcnt", "-mlzcnt"),
    ("movbe", "-mmovbe"),
    ("f16c", "-mf16c"),
];

/// Return the compiler flags for the target features, in a fixed order.
/// `features`: the names of the Rust target features, as listed in `CARGO_CFG_TARGET_FEATURE`.
/// Features that have no entry in the table are ignored.
pub fn compiler_flags<'a>(features: impl IntoIterator<Item = &'a str>) -> Vec<&'static str> {
    let features: Vec<&str> = features.into_iter().collect();
    FEATURE_FLAGS
        .iter()
        .filter(|(feature, _)| features.contains(feature))
        .map(|&(_, flag)| flag)
        .collect()
}
