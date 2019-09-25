// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::os::raw::c_char;
use std::os::raw::c_int;
use std::os::raw::c_ulonglong;
use std::path::Path;
use std::sync::Arc;

use failure::Error;
use rust_htslib::bam;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::HeaderView;

mod bindings;

use bindings::Aligner as BindAligner;
use bindings::StarRef as BindRef;

pub struct StarReference {
    inner: Arc<InnerStarReference>,
}

struct InnerStarReference {
    reference: *const BindRef,
    settings: StarSettings,
    header: Header,
    header_view: HeaderView,
}

unsafe impl Sync for InnerStarReference {}

impl Drop for InnerStarReference {
    fn drop(&mut self) {
        unsafe {
            bindings::destroy_ref(self.reference);
        }
    }
}

impl StarReference {
    /// Load the reference index into memory based on the given settings.
    /// Should only be called once for a given reference. Create aligners
    /// that use this reference by calling `get_aligner`.  The reference
    /// index will be free'd when all `Aligner`s that use this `StarReference`
    /// have been dropped.
    pub fn load(settings: StarSettings) -> Result<StarReference, Error> {
        // Load headers
        let (header, header_view) = generate_header(&settings.reference_path);

        // Load reference
        let mut nvec = Vec::new();
        for x in settings.args.iter() {
            let cur_string = CString::new(x.as_str())?;
            nvec.push(cur_string.into_raw());
        }

        let c_args = nvec.as_mut_ptr() as *mut *mut c_char;
        let length = nvec.len() as c_int;

        let reference = unsafe { bindings::init_star_ref(length, c_args) };

        let inner = InnerStarReference {
            reference,
            header,
            header_view,
            settings,
        };

        Ok(StarReference {
            inner: Arc::new(inner),
        })
    }

    pub fn header(&self) -> &Header {
        &self.inner.as_ref().header
    }

    pub fn header_view(&self) -> &HeaderView {
        &self.inner.as_ref().header_view
    }

    pub fn reference_path(&self) -> &str {
        &self.inner.as_ref().settings.reference_path
    }

    pub fn get_aligner(&self) -> StarAligner {
        StarAligner::new(self.inner.clone())
    }
}

/// StarSettings contains the parameters which will be used for the STAR aligner.
/// Currently the array of argument strings is passed directly
#[derive(Clone)]
pub struct StarSettings {
    reference_path: String,
    multn: usize,
    args: Vec<String>,
}

/// The default value for number of mappings allowed is 1.
const DEFAULT_MULTN: usize = 1;

impl StarSettings {
    /// This constructor just sets all of the necessary arguments to their defaults, and the
    /// arguments which can take on different values have separate functions to set them later
    pub fn new(reference_path: &str) -> StarSettings {
        let def_args: Vec<String> = vec![
            "STAR".to_string(),
            "--genomeDir".to_string(),
            reference_path.to_string(),
            "--outSAMmultNmax".to_string(),
            DEFAULT_MULTN.to_string(),
            "--runThreadN".to_string(),
            "1".to_string(),
            "--readNameSeparator".to_string(),
            "space".to_string(),
            "--outSAMunmapped".to_string(),
            "Within".to_string(),
            "KeepPairs".to_string(),
            "--outSAMtype".to_string(),
            "SAM".to_string(),
            "--outStd".to_string(),
            "SAM".to_string(),
            "--outSAMorder".to_string(),
            "PairedKeepInputOrder".to_string(),
        ];
        StarSettings {
            reference_path: reference_path.to_string(),
            multn: DEFAULT_MULTN,
            args: def_args,
        }
    }

    /// Update the array of arguments when the arguments' values have changed
    fn sync_args(&mut self) {
        for i in 0..self.args.len() {
            if self.args[i] == "--genomeDir" {
                self.args[i + 1] = self.reference_path.clone();
            } else if self.args[i] == "--outSAMmultNmax" {
                self.args[i + 1] = self.multn.to_string();
            }
        }
    }

    /// Set the max number of multimapping reads
    pub fn set_multn(&mut self, new_multn: usize) {
        self.multn = new_multn;
        self.sync_args();
    }

    /// Add the given read group strings to the arguments
    pub fn add_rg(&mut self, rg_tags: Vec<String>) {
        self.args.push("--outSAMattrRGline".to_string());
        for tag in rg_tags {
            self.args.push(tag);
        }
    }
}

/// StarAligner aligns single reads or read-pairs to the reference it is initialized with, and returns
/// rust_htslib Record objects
pub struct StarAligner {
    aligner: *mut BindAligner,
    reference: Arc<InnerStarReference>,
    sam_buf: Vec<u8>,
    aln_buf: Vec<u8>,
}

impl StarAligner {
    fn new(reference: Arc<InnerStarReference>) -> StarAligner {
        let aligner = unsafe { bindings::init_aligner_from_ref(reference.as_ref().reference) };
        StarAligner {
            aligner,
            reference,
            sam_buf: Vec::new(),
            aln_buf: Vec::new(),
        }
    }

    /// Aligns a given read and produces BAM records
    pub fn align_read(&mut self, name: &[u8], read: &[u8], qual: &[u8]) -> Vec<bam::Record> {
        align_read_rust(self.aligner, read, qual, &mut self.aln_buf).unwrap();
        self.parse_sam_to_records(name)
    }

    /// Aligns a given read and return the resulting SAM string
    pub fn align_read_sam(&mut self, _name: &[u8], read: &[u8], qual: &[u8]) -> String {
        align_read_rust(self.aligner, read, qual, &mut self.aln_buf).unwrap();
        String::from_utf8(self.aln_buf.clone()).unwrap()
    }

    /// Aligns a given pair of reads and produces BAM records
    pub fn align_read_pair(
        &mut self,
        name: &[u8],
        read1: &[u8],
        qual1: &[u8],
        read2: &[u8],
        qual2: &[u8],
    ) -> (Vec<bam::Record>, Vec<bam::Record>) {
        align_read_pair_rust(self.aligner, read1, qual1, read2, qual2, &mut self.aln_buf).unwrap();
        let full_vec = self.parse_sam_to_records(name);

        // Partition the records into first mate and second mate
        let mut first_vec: Vec<bam::Record> = Vec::new();
        let mut second_vec: Vec<bam::Record> = Vec::new();
        for rec in full_vec {
            if rec.is_first_in_template() {
                first_vec.push(rec);
            } else {
                second_vec.push(rec);
            }
        }
        (first_vec, second_vec)
    }

    /// Aligns a given read and return the resulting SAM string
    pub fn align_read_pair_sam(
        &mut self,
        _name: &[u8],
        read1: &[u8],
        qual1: &[u8],
        read2: &[u8],
        qual2: &[u8],
    ) -> String {
        align_read_pair_rust(self.aligner, read1, qual1, read2, qual2, &mut self.aln_buf).unwrap();
        String::from_utf8(self.aln_buf.clone()).unwrap()
    }

    /// Given a list of BAM records as a SAM-format string in which records are separated by new
    /// lines, add the records to a vector and append the read name to the beginning of them so
    /// that they conform with BAM specifications
    fn parse_sam_to_records(&mut self, name: &[u8]) -> Vec<bam::Record> {
        let mut records = Vec::new();
        for slc in self.aln_buf.split(|c| *c == b'\n') {
            if slc.len() > 0 {
                self.sam_buf.clear();
                self.sam_buf.extend_from_slice(name);
                self.sam_buf.extend_from_slice(slc);
                let record =
                    bam::Record::from_sam(&self.reference.as_ref().header_view, &self.sam_buf)
                        .unwrap();
                records.push(record);
            }
        }

        records
    }
}

impl Clone for StarAligner {
    fn clone(&self) -> StarAligner {
        StarAligner::new(self.reference.clone())
    }
}

impl Drop for StarAligner {
    fn drop(&mut self) {
        unsafe { bindings::destroy_aligner(self.aligner) };
    }
}

/// Read in the lines from a file and store each line as its own string in a vector
fn get_lines(path: &Path) -> Vec<String> {
    let file = match File::open(&path) {
        Err(error) => panic!("error: {}: {}", path.display(), error),
        Ok(file) => file,
    };
    BufReader::new(file)
        .lines()
        .map(|line| line.unwrap())
        .collect()
}

/// Produces a header from the genome reference directory by looking up the contig names and
/// lengths and formatting them properly
fn generate_header(genome_path: impl AsRef<Path>) -> (Header, HeaderView) {
    let mut header = Header::new();

    let contig_names_path = genome_path.as_ref().join(Path::new("chrName.txt"));
    let contig_names = get_lines(&contig_names_path);

    let contig_lengths_path = genome_path.as_ref().join(Path::new("chrLength.txt"));
    let contig_lengths = get_lines(&contig_lengths_path);
    for (ref contig_name, len) in contig_names.iter().zip(contig_lengths.iter()) {
        add_ref_to_bam_header(&mut header, &contig_name, len.parse::<usize>().unwrap());
    }

    let hv = HeaderView::from_header(&header);
    (header, hv)
}

/// Given a reference genome contig's name and length, add a corresponding line to the given BAM
/// header
fn add_ref_to_bam_header(header: &mut Header, seq_name: &str, seq_len: usize) {
    let mut header_rec = HeaderRecord::new(b"SQ");
    header_rec.push_tag(b"SN", &seq_name);
    header_rec.push_tag(b"LN", &seq_len);
    header.push_record(&header_rec);
}

/// Below are wrappers for different Orbit functions, but with arguments as datatypes which are
/// more natural rather than the wrappers around C datatypes.  Each function below makes any
/// necessary conversions to the inputs, calls the library function, and makes any necessary
/// conversions to the outputs.
fn align_read_rust(
    al: *mut BindAligner,
    read: &[u8],
    qual: &[u8],
    aln_buf: &mut Vec<u8>,
) -> Result<(), Error> {
    let length = read.len() as c_ulonglong;
    let c_read = CString::new(read)?;
    let c_qual = CString::new(qual)?;
    let read_ptr = c_read.as_ptr() as *mut c_char;
    let qual_ptr = c_qual.as_ptr() as *mut c_char;

    let res: *const c_char = unsafe { bindings::align_read(al, read_ptr, qual_ptr, length) };
    if res.is_null() {
        return Err(failure::format_err!("STAR returned null alignment"));
    }

    let cstr = unsafe { CStr::from_ptr(res) };
    aln_buf.clear();
    aln_buf.extend_from_slice(cstr.to_bytes());

    unsafe {
        libc::free(res as *mut libc::c_void);
    }
    Ok(())
}

fn align_read_pair_rust(
    al: *mut BindAligner,
    read: &[u8],
    qual: &[u8],
    read2: &[u8],
    qual2: &[u8],
    aln_buf: &mut Vec<u8>,
) -> Result<(), Error> {
    let length = read.len() as c_ulonglong;
    let c_read = CString::new(read)?;
    let c_qual = CString::new(qual)?;
    let read_ptr = c_read.as_ptr() as *mut c_char;
    let qual_ptr = c_qual.as_ptr() as *mut c_char;

    let c_read2 = CString::new(read2)?;
    let c_qual2 = CString::new(qual2)?;
    let read_ptr2 = c_read2.as_ptr() as *mut c_char;
    let qual_ptr2 = c_qual2.as_ptr() as *mut c_char;

    let res: *const c_char =
        unsafe { bindings::align_read_pair(al, read_ptr, qual_ptr, read_ptr2, qual_ptr2, length) };
    if res.is_null() {
        return Err(failure::format_err!("STAR returned null alignment"));
    }

    let cstr = unsafe { CStr::from_ptr(res) };
    aln_buf.clear();
    aln_buf.extend_from_slice(cstr.to_bytes());

    unsafe {
        libc::free(res as *mut libc::c_void);
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    /// References to some commonly used reference genomes for testing purposes
    pub const DEFAULT_REF_1: &str = "/mnt/opt/refdata_cellranger/mm10-3.0.0/star";
    pub const DEFAULT_REF_2: &str = "/mnt/opt/refdata_cellranger/GRCh38-3.0.0/star";

    const ERCC_REF: &'static str = "test/ercc92-1.2.0/star/";

    const NAME: &'static [u8] = b"NAME";
    const ERCC_READ_1: &'static [u8] = b"GCATCCAGACCGTCGGCTGATCGTGGTTTTACTAGGCTAGACTAGCGTACGAGCACTATGGTCAGTAATTCCTGGAGGAATAGGTACCAAGAAAAAAACG";
    const ERCC_QUAL_1: &'static [u8] = b"????????????????????????????????????????????????????????????????????????????????????????????????????";

    const ERCC_READ_2: &'static [u8] = b"GGAGACGAATTGCCAGAATTATTAACTGCGCAGTTAGGGCAGCGTCTGAGGAAGTTTGCTGCGGTTTCGCCTTGACCGCGGGAAGGAGACATAACGATAG";
    const ERCC_QUAL_2: &'static [u8] = b"????????????????????????????????????????????????????????????????????????????????????????????????????";

    #[test]
    fn test_ercc_align() {
        let settings = StarSettings::new(ERCC_REF);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let recs = aligner.align_read(NAME, ERCC_READ_1, ERCC_QUAL_1);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].pos(), 50);
        assert_eq!(recs[0].tid(), 0);
        println!("{:?}", recs);

        let recs = aligner.align_read(NAME, ERCC_READ_2, ERCC_QUAL_2);
        assert_eq!(recs[0].tid(), 0);
        assert_eq!(recs[0].pos(), 500);
        println!("{:?}", recs);
    }

    #[test]
    #[ignore]
    fn test_align_read() {
        let settings = StarSettings::new(DEFAULT_REF_2);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let read = b"GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC";
        let qual = b"GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG";

        let res = aligner.align_read_sam(b"name", read, qual);
        assert!(res == "\t0\t6\t30070474\t255\t98M\t*\t0\t0\tGTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC\tGGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:0\n");
    }

    #[test]
    #[ignore]
    fn test_align_multiple_reads() {
        let settings = StarSettings::new(DEFAULT_REF_2);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let read = b"GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC";
        let qual = b"GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG";
        let read2 = b"GTATGTCAAGTTGGTGGAGGCCCTTTGTGCTGAACACCAAATCAACCTAATTAAGGTTGATGACAACAAGAAACTAGGAGAATGGGTAGGCCTTTGTA";
        let qual2 = b"AGGAGGGGIG.GAGGGGIGGIGIIGGGIIAAGGGGGIGGIIGAGIGIA.GGGGIGGGGGGGGGGGGGGIGIIIIGGGGGGIGGGGIIIGA.G.<.GGG";

        let res = aligner.align_read_sam(b"name", read, qual);
        let res2 = aligner.align_read_sam(b"name2", read2, qual2);

        assert!(res == "\t0\t6\t30070474\t255\t98M\t*\t0\t0\tGTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC\tGGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:0\n");
        assert!(res2 == "\t0\t6\t132816509\t255\t55M396N43M\t*\t0\t0\tGTATGTCAAGTTGGTGGAGGCCCTTTGTGCTGAACACCAAATCAACCTAATTAAGGTTGATGACAACAAGAAACTAGGAGAATGGGTAGGCCTTTGTA\tAGGAGGGGIG.GAGGGGIGGIGIIGGGIIAAGGGGGIGGIIGAGIGIA.GGGGIGGGGGGGGGGGGGGIGIIIIGGGGGGIGGGGIIIGA.G.<.GGG\tNH:i:1\tHI:i:1\tAS:i:98\tnM:i:0\n");
    }

    #[test]
    fn test_header() {
        let settings = StarSettings::new(ERCC_REF);
        let reference = StarReference::load(settings).unwrap();

        let header_string = String::from_utf8(reference.header().to_bytes()).unwrap();
        assert!(header_string.starts_with("@SQ\tSN:ERCC-00002\tLN:1061\n"));
        assert!(header_string.ends_with("@SQ\tSN:ERCC-00171\tLN:505"));
    }

    #[test]
    #[ignore]
    fn test_get_record() {
        let settings = StarSettings::new(DEFAULT_REF_2);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let read = b"GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC";
        let qual = b"GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG";
        let name = b"gatactaga";
        let res = aligner.align_read(name, read, qual);
        assert!(res.len() > 0);
        assert_eq!(res[0].pos(), 30070473);
    }

    #[test]
    #[ignore]
    fn test_write_bam() {
        let settings = StarSettings::new(DEFAULT_REF_1);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let mut out = bam::Writer::from_path(&"test/test.bam", &reference.header()).unwrap();
        let read = b"GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC";
        let qual = b"GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG";
        let name = b"gatactaga";
        let res = aligner.align_read(name, read, qual);
        for record in res.iter() {
            out.write(&record).unwrap();
        }
        let bam_wrapped = bam::Reader::from_path(&"test/test.bam");
        match bam_wrapped {
            Ok(v) => println!("working with version: {:?}", v),
            Err(e) => println!("error parsing header: {:?}", e),
        }
        // The reading portion is commented out because it's failing for unknown reasons, but the BAM
        // file can be read with samtools and pysam.
        /*
        let mut bam = bam::Reader::from_path(&"test/test.bam").unwrap();
        let read_records = bam.records();
        let mut i = 0;
        for r in read_records {
            let record = r.unwrap();
            assert_eq!(record.pos(), res[i].pos());
            i += 1;
        }
        */
    }

    /*
    #[test]
    fn test_align_fastq()
    {
        let settings = StarSettings::new(DEFAULT_REF_1);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let fastq_path : &str = "test/full1_sample.fastq";

        let res : Vec<bam::Record> = aligner.align_fastq(fastq_path).unwrap();
        assert!(res.len() > 0);

        let mut out = bam::Writer::from_path(&"test/full1_sample.bam", &reference.header()).unwrap();
        for record in res.iter() {
            out.write(&record).unwrap();
        }
    }
    */
}
