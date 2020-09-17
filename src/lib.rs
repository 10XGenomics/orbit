// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::os::raw::c_char;
use std::os::raw::c_int;
use std::path::Path;
use std::sync::Arc;

use failure::Error;
use rust_htslib::bam;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::HeaderView;
use star_sys::{self as bindings, Aligner as BindAligner, StarRef as BindRef};

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

        // recover stray CStrings to prevent leaked memory
        nvec.into_iter().for_each(|ptr| unsafe {
            CString::from_raw(ptr);
        });

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
    args: Vec<String>,
}

impl StarSettings {
    /// This constructor just sets all of the necessary arguments to their defaults, and the
    /// arguments which can take on different values have separate functions to set them later
    pub fn new(reference_path: &str) -> StarSettings {
        let def_args: Vec<String> = vec![
            "STAR".to_string(),
            "--genomeDir".to_string(),
            reference_path.to_string(),
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
            args: def_args,
        }
    }

    /// Add an argument to pass to STAR.
    pub fn arg(mut self, arg: &str) -> Self {
        self.args.push(arg.to_string());
        self
    }

    /// Add the given read group strings to the arguments
    pub fn add_rg(&mut self, rg_tags: Vec<String>) {
        self.args.push("--outSAMattrRGline".to_string());
        for tag in rg_tags {
            self.args.push(tag);
        }
    }

    /// Provide an estimate of the memory usage (in bytes) of the reference
    pub fn est_mem(&self) -> Result<usize, Error> {
        let refpath = Path::new(&self.reference_path);
        ["Genome", "SA", "SAindex"]
            .iter()
            .fold(Ok(0usize), |acc, file| {
                Ok(acc? + std::fs::metadata(refpath.join(file))?.len() as usize)
            })
    }
}

/// StarAligner aligns single reads or read-pairs to the reference it is initialized with, and returns
/// rust_htslib Record objects
pub struct StarAligner {
    aligner: *mut BindAligner,
    reference: Arc<InnerStarReference>,
    sam_buf: Vec<u8>,
    aln_buf: Vec<u8>,
    fastq1: Vec<u8>,
    fastq2: Vec<u8>,
    header_view: HeaderView,
}

unsafe impl Send for StarAligner {}

enum AlignedRecords<'a> {
    Read1(&'a mut Vec<bam::Record>),
    Read2(&'a mut Vec<bam::Record>),
}

impl StarAligner {
    fn new(reference: Arc<InnerStarReference>) -> StarAligner {
        let aligner = unsafe { bindings::init_aligner_from_ref(reference.as_ref().reference) };
        let header_view = reference.as_ref().header_view.clone();

        StarAligner {
            aligner,
            reference,
            sam_buf: Vec::new(),
            aln_buf: Vec::new(),
            fastq1: Vec::new(),
            fastq2: Vec::new(),
            header_view,
        }
    }

    /// Prepare a read and qual for processing by formatting as FASTQ
    fn prepare_fastq(buf: &mut Vec<u8>, _name: &[u8], read: &[u8], qual: &[u8]) {
        buf.clear();
        // for now, we fixup the read name on the backend, avoid the copy for now
        buf.extend_from_slice(b"@a\n");
        buf.extend_from_slice(read);
        buf.extend_from_slice(b"\n+\n");
        buf.extend_from_slice(qual);
        buf.push(b'\n');
        buf.push(b'\0');
    }

    fn empty_record(name: &[u8], read: &[u8], qual: &[u8]) -> bam::Record {
        let mut rec = bam::Record::new();
        rec.set_tid(-1);
        rec.set_pos(-1);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        rec.set(name, None, read, qual);
        rec.set_unmapped();
        rec
    }

    fn empty_records<'a>(
        name: &[u8],
        read: &[u8],
        qual: &[u8],
        alns: AlignedRecords<'a>,
    ) -> Vec<bam::Record> {
        use AlignedRecords::*;
        let (alns, aln_is_first_in_template) = match alns {
            Read1(alns) => (alns, true),
            Read2(alns) => (alns, false),
        };
        let mut recs = vec![];
        recs.resize_with(alns.len(), || Self::empty_record(name, read, qual));
        for (aln, rec) in alns.iter_mut().zip(&mut recs) {
            aln.set_paired();
            rec.set_paired();
            if aln.is_secondary() {
                rec.set_secondary();
            }
            if aln.is_supplementary() {
                rec.set_supplementary();
            }
            if aln_is_first_in_template {
                // aln is first in template and rec is an empty unmapped record
                aln.set_first_in_template();
                rec.set_last_in_template();
                aln.set_mate_unmapped();
            } else {
                // rec is first in template (and an empty unmapped record),
                //   and aln is last in template
                rec.set_first_in_template();
                aln.set_last_in_template();
                if aln.is_unmapped() {
                    rec.set_mate_unmapped();
                }
            }
        }
        recs
    }

    /// Aligns a given read and produces BAM records
    pub fn align_read(&mut self, name: &[u8], read: &[u8], qual: &[u8]) -> Vec<bam::Record> {
        // STAR will throw an error on empty reads - so just construct an empty record.
        if read.len() == 0 {
            // Make an unmapped record and return it
            return vec![Self::empty_record(name, read, qual)];
        }

        Self::prepare_fastq(&mut self.fastq1, name, read, qual);
        align_read_rust(self.aligner, self.fastq1.as_slice(), &mut self.aln_buf).unwrap();
        self.parse_sam_to_records(name)
    }

    /// Aligns a given read and return the resulting SAM string
    pub fn align_read_sam(&mut self, name: &[u8], read: &[u8], qual: &[u8]) -> String {
        Self::prepare_fastq(&mut self.fastq1, name, read, qual);
        align_read_rust(self.aligner, self.fastq1.as_slice(), &mut self.aln_buf).unwrap();
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
        use AlignedRecords::*;
        if read1.len() == 0 {
            let mut recs2 = self.align_read(name, read2, qual2);
            let recs1 = Self::empty_records(name, read1, qual1, Read2(&mut recs2));
            return (recs1, recs2);
        } else if read2.len() == 0 {
            let mut recs1 = self.align_read(name, read1, qual1);
            let recs2 = Self::empty_records(name, read2, qual2, Read1(&mut recs1));
            return (recs1, recs2);
        }
        Self::prepare_fastq(&mut self.fastq1, name, read1, qual1);
        Self::prepare_fastq(&mut self.fastq2, name, read2, qual2);
        align_read_pair_rust(
            self.aligner,
            self.fastq1.as_slice(),
            self.fastq2.as_slice(),
            &mut self.aln_buf,
        )
        .unwrap();
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
        name: &[u8],
        read1: &[u8],
        qual1: &[u8],
        read2: &[u8],
        qual2: &[u8],
    ) -> String {
        Self::prepare_fastq(&mut self.fastq1, name, read1, qual1);
        Self::prepare_fastq(&mut self.fastq2, name, read2, qual2);
        align_read_pair_rust(
            self.aligner,
            self.fastq1.as_slice(),
            self.fastq2.as_slice(),
            &mut self.aln_buf,
        )
        .unwrap();
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
                let record = bam::Record::from_sam(&mut self.header_view, &self.sam_buf).unwrap();
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
fn align_read_rust(al: *mut BindAligner, fastq: &[u8], aln_buf: &mut Vec<u8>) -> Result<(), Error> {
    let fastq = CStr::from_bytes_with_nul(fastq)?;
    let res: *const c_char = unsafe { bindings::align_read(al, fastq.as_ptr()) };
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
    fastq1: &[u8],
    fastq2: &[u8],
    aln_buf: &mut Vec<u8>,
) -> Result<(), Error> {
    let fastq1 = CStr::from_bytes_with_nul(fastq1)?;
    let fastq2 = CStr::from_bytes_with_nul(fastq2)?;
    let res: *const c_char =
        unsafe { bindings::align_read_pair(al, fastq1.as_ptr(), fastq2.as_ptr()) };
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
    use rust_htslib::bam::Read;

    /// References to some commonly used reference genomes for testing purposes
    pub const DEFAULT_REF_1: &str = "/mnt/opt/refdata_cellranger/mm10-3.0.0/star";
    pub const DEFAULT_REF_2: &str = "/mnt/opt/refdata_cellranger/GRCh38-3.0.0/star";
    pub const DEFAULT_REF_3: &str = "/mnt/opt/refdata_cellranger/GRCh38-1.2.0/star";
    pub const DEFAULT_REF_4: &str = "//mnt/opt/refdata_cellranger/hg19_and_mm10-3.0.0/star";

    const ERCC_REF: &'static str = "test/ercc92-1.2.0/star/";

    const NAME: &'static [u8] = b"NAME";
    const ERCC_READ_1: &'static [u8] = b"GCATCCAGACCGTCGGCTGATCGTGGTTTTACTAGGCTAGACTAGCGTACGAGCACTATGGTCAGTAATTCCTGGAGGAATAGGTACCAAGAAAAAAACG";
    const ERCC_QUAL_1: &'static [u8] = b"????????????????????????????????????????????????????????????????????????????????????????????????????";

    const ERCC_READ_2: &'static [u8] = b"GGAGACGAATTGCCAGAATTATTAACTGCGCAGTTAGGGCAGCGTCTGAGGAAGTTTGCTGCGGTTTCGCCTTGACCGCGGGAAGGAGACATAACGATAG";
    const ERCC_QUAL_2: &'static [u8] = b"????????????????????????????????????????????????????????????????????????????????????????????????????";

    const ERCC_READ_3: &'static [u8] = b"AACTTAATGGACGGG";
    const ERCC_QUAL_3: &'static [u8] = b"???????????????";

    const ERCC_READ_4: &'static [u8] = b"AATCCACTCAATAAATCTAAAAAC";
    const ERCC_QUAL_4: &'static [u8] = b"????????????????????????";

    fn have_refs() -> bool {
        Path::new("/mnt/opt/refdata_cellranger").exists()
    }

    #[test]
    fn test_empty_tiny_reads() {
        let settings = StarSettings::new(ERCC_REF);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let recs = aligner.align_read(b"a", b"", b"");
        println!("{:?}", recs);
        let recs = aligner.align_read(b"b", b"A", b"?");
        println!("{:?}", recs);
        let (recs1, recs2) = aligner.align_read_pair(b"a", b"", b"", b"", b"");
        println!("{:?}, {:?}", recs1, recs2);
        let (recs1, recs2) = aligner.align_read_pair(b"b", b"A", b"?", b"", b"");
        println!("{:?}, {:?}", recs1, recs2);
        let (recs1, recs2) = aligner.align_read_pair(b"c", b"", b"", b"C", b"?");
        println!("{:?}, {:?}", recs1, recs2);
        let (recs1, recs2) = aligner.align_read_pair(b"d", b"A", b"?", b"C", b"?");
        println!("{:?}, {:?}", recs1, recs2);
    }

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
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].tid(), 0);
        assert_eq!(recs[0].pos(), 500);
        println!("{:?}", recs);

        let recs = aligner.align_read(NAME, ERCC_READ_3, ERCC_QUAL_3);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].flags(), 0);
        assert_eq!(recs[0].tid(), 39);
        assert_eq!(recs[0].pos(), 27);
        assert_eq!(recs[0].mapq(), 3);
        assert_eq!(recs[1].flags(), 0x110); // REVERSE,SECONDARY
        assert_eq!(recs[1].tid(), 72);
        assert_eq!(recs[1].pos(), 553);
        assert_eq!(recs[1].mapq(), 3);
        println!("{:?}", recs);

        let recs = aligner.align_read(NAME, ERCC_READ_4, ERCC_QUAL_4);
        println!("{:?}", recs);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].flags(), 0);
        assert_eq!(recs[0].tid(), 72);
        assert_eq!(recs[0].pos(), 492);
        assert_eq!(recs[0].mapq(), 3);
        assert_eq!(recs[1].flags(), 0x100); // SECONDARY
        assert_eq!(recs[1].tid(), 72);
        assert_eq!(recs[1].pos(), 607);
        assert_eq!(recs[1].mapq(), 3);
    }

    #[test]
    fn test_transcriptome_min_score() {
        let settings = StarSettings::new(ERCC_REF).arg("--outFilterScoreMin=20");
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let recs = aligner.align_read(NAME, ERCC_READ_3, ERCC_QUAL_3);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].flags(), 4); // UNMAP
        assert_eq!(recs[0].tid(), -1);
        assert_eq!(recs[0].pos(), -1);
        assert_eq!(recs[0].mapq(), 0);
        println!("{:?}", recs);
    }

    #[test]
    fn test_ercc_align_zero_len() {
        let settings = StarSettings::new(ERCC_REF);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let zero_base_recs = aligner.align_read(NAME, b"", b"");
        assert_eq!(zero_base_recs.len(), 1);
        assert!(zero_base_recs[0].is_unmapped());

        println!(
            "zero: {:?}, flags {:#b}",
            zero_base_recs,
            zero_base_recs[0].flags()
        );

        let one_base_recs = aligner.align_read(NAME, b"G", b"D");
        assert_eq!(one_base_recs.len(), 1);
        assert!(one_base_recs[0].is_unmapped());

        println!(
            "one: {:?}, flags {:#b}",
            one_base_recs,
            one_base_recs[0].flags()
        );

        assert_eq!(zero_base_recs[0].tid(), one_base_recs[0].tid());
        assert_eq!(zero_base_recs[0].pos(), one_base_recs[0].pos());
        assert_eq!(zero_base_recs[0].mpos(), one_base_recs[0].mpos());
        assert_eq!(zero_base_recs[0].mtid(), one_base_recs[0].mtid());
        assert_eq!(zero_base_recs[0].flags(), one_base_recs[0].flags());
    }

    #[test]
    fn test_multithreaded_alignment() {
        let settings = StarSettings::new(ERCC_REF);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner1 = reference.get_aligner();
        let mut aligner2 = reference.get_aligner();

        let t1 = std::thread::spawn(move || {
            for _ in 0..100000 {
                let recs = aligner1.align_read(NAME, ERCC_READ_1, ERCC_QUAL_1);
                assert_eq!(recs.len(), 1);
                assert_eq!(recs[0].pos(), 50);
                assert_eq!(recs[0].tid(), 0);
            }
        });

        let t2 = std::thread::spawn(move || {
            for _ in 0..100000 {
                let recs = aligner2.align_read(NAME, ERCC_READ_2, ERCC_QUAL_2);
                assert_eq!(recs.len(), 1);
                assert_eq!(recs[0].pos(), 500);
                assert_eq!(recs[0].tid(), 0);
            }
        });

        assert!(t1.join().is_ok());
        assert!(t2.join().is_ok());
    }

    #[test]
    fn test_align_read() {
        if !have_refs() {
            return;
        }

        let settings = StarSettings::new(DEFAULT_REF_2);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let read = b"GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC";
        let qual = b"GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG";

        let res = aligner.align_read_sam(b"name", read, qual);
        assert!(res == "\t0\t6\t30070474\t255\t98M\t*\t0\t0\tGTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC\tGGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:0\n");
    }

    // random selection of genome sequences
    const SEQ: &[&[u8]] = &[
        b"TTCCGTATATTAAACCAACCCTGTATGCTTGCGATGAAGCCTATTTGATCCTGTTGGATGATTGTTTTAATGTGCTCTTGGATTCGGTTTTCCAGAATTTACTGAGTATTTTTGCATCGA",
        b"GTAGTGAGAGAGGCAAACAGTCGGGATTGGATACAAGGAAGTCACTCTGATTACCAGCCTTTGAAGGTCAACAGGGGGACGGGTGAGAAGGAATGGGAGGCCCAGCGGCCTGGAGTCCTG",
        b"CAAAAGAAGCCAACAGAAAAAACTGAAAATTGAGAAATCTTTTTATCATGATTGAATTGTTCTAAAATTGGACTGTGATGGAGACACATCTTTACACTTAGACACACAAAAAAGCAAAGA",
        b"CGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATGATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTT",
        b"CTTGATGAGGTCACCCGGCATGAAAGCTTTCACCCCGTTGGAGGGCTTCCCTAGGGAGTCCTGGGTTATTCCGTAGACATGCGACCCCTCTACTGCTCCTTCCCTCATGGAGCCTCCCCT",
        b"GCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAGGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCAGAGGTGA",
        b"CTCCAGGGCGTCGCCCAGGTTGAAGTCGTCACTCGCTGAGGGGGGCGGGGACAAGATTATCATAAATATGCAAATTTAGCGCCTCATTTACATATACAGCATTTACATCCAATCCGGGCC",
        b"CCGCCCCCAGGCCCCACAGCTCGAGCTGATTTGAATATGCAGATGAGCAGACTCCGCCCCCGGGCACTGCCTTTAACCCCGCCCATGAGCCTGACGCTGAGTTCCGGCTTCCAAGCCCCG",
        b"GCCCAGAGACAAGGACTTGCAGGGCAACAAAGCAGCATCATCCCACTATTCCAGGGGAGGTGCTAAATACGAGGGTGAGGCTGTCAAGCAGTCCCTGGTG",
        b"GGCCCTTCCTCCCAGCCCAGACTCCTACATCCCAAACTTGAGCCATGGCACACATGCTGGGCACTTACTCTGTGCATAGCAGAGGGAGCTGAGCTGCATC",
        b"TTCTAGCCCCTCCAACAGTTCTTCATGTAGCACTGCTTTTGGACAGTGACCCATGATATTGCCTTCGTTTGCCCCTCTGACTGTTGGAATCCCACGCACC",
        b"GATGAGAAGCTCTGGCAGCTGGTAGCCATGGCGAAGATAGAGAGGTTCTCGTATGGGCAGCTGATCTCAAAAGATTTTGGAGAGTCACCCTTCATCATGT",
        b"TGTAATTTTCTAAGTGATAAGAGGACTAGGAGCATCTTTTGTTCTAATGAGGTGACTCTGGGTGAACTCCTGGCTGGGTCCTGGATGGGGGCTGGTCAAG"
    ];

    // random selection of references
    const REF_SET: &[&str] = &[
        "hg19_and_mm10-3.0.0",
        "Rattus_norvegicus.Rnor_6.0",
        "Rattus_norvegicus.Rnor_6.0_EGFP",
        "Danio_rerio_GRCz11.99_eGFP_DsRed",
        "ASM130575v1-1.0.0",
        "EquCab2-1.0.0",
        "mm10-3.0.0",
        "Mmul_8.0.1",
        "Sscrofa11.1-1.0.0",
        "donkey_3611_h1_MAKER",
        "ppatens3_3-1.0.0",
        "bdgp6-1.0.0",
    ];

    #[test]
    fn test_scan_references() {
        if !have_refs() {
            return;
        }

        let path: String = "/mnt/opt/refdata_cellranger".to_string();

        for ref_test in REF_SET {
            let p = path.clone() + "/" + ref_test + "/star";

            let settings = StarSettings::new(&p);
            let reference = StarReference::load(settings).unwrap();
            let mut aligner = reference.get_aligner();

            for read_test in SEQ {
                let qual = vec![b'I'; read_test.len()];
                let res = aligner.align_read_sam(b"asdf", read_test, &qual);
                println!("res: {}", res);
            }
        }
    }

    #[test]
    fn test_align_read_hgmm() {
        if !have_refs() {
            return;
        }

        let settings = StarSettings::new(DEFAULT_REF_4);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let read = b"GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC";
        let qual = b"GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG";

        let res = aligner.align_read_sam(b"name", read, qual);
        println!("res: {}", res);
        assert!(res == "\t0\thg19_6\t30038251\t255\t98M\t*\t0\t0\tGTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC\tGGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:0\n");
    }

    #[test]
    fn test_align_read_pair() {
        if !have_refs() {
            return;
        }

        let settings = StarSettings::new(DEFAULT_REF_3);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let read1 = b"GTATAAACAAAAATCCTGTCTCTAAAATGTAACTGATGACTAGCTGACATGCAGCTACAAAGACCTTAGTTCCTTTAAAAACAATATCCAATCAAATCAGAATTGCCCCTCATGCTTCACATAA";
        let qual1 = b"-FFF8FFFF8FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8FF8FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        let read2 = b"GAGTTGTGTAAAGTGGCCAAACATCAACAACAACAACAAAAAACACAAACAGGAAAGAGCAATTGGGTAAAGCTGTACACTGCTCTTTTAAAATACTATTATGAGATGGACATTTATGTGAAGCATGAGGGGCAATTCTGATTTGATTGG";
        let qual2 = b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFF8FFFFFF-FFFFFFF8";

        let (recs1, recs2) = aligner.align_read_pair(b"name", read1, qual1, read2, qual2);

        /*
        163	1	24647211	255	150M	=	24647324	237	GAGTTGTGTAAAGTGGCCAAACATCAACAACAACAACAAAAAACACAAACAGGAAAGAGCAATTGGGTAAAGCTGTACACTGCTCTTTTAAAATACTATTATGAGATGGACATTTATGTGAAGCATGAGGGGCAATTCTGATTTGATTGG	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFF8FFFFFF-FFFFFFF8	NH:i:1	HI:i:1	AS:i:272	nM:i:0
        83	1	24647324	255	124M	=	24647211	-237	TTATGTGAAGCATGAGGGGCAATTCTGATTTGATTGGATATTGTTTTTAAAGGAACTAAGGTCTTTGTAGCTGCATGTCAGCTAGTCATCAGTTACATTTTAGAGACAGGATTTTTGTTTATAC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8FF8FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8FFFF8FFF-	NH:i:1	HI:i:1	AS:i:272	nM:i:0
        */

        assert_eq!(recs1.len(), 1);
        assert_eq!(recs1[0].flags(), 83);
        assert_eq!(recs1[0].tid(), 0);
        assert_eq!(recs1[0].pos(), 24647323);

        assert_eq!(recs2.len(), 1);
        assert_eq!(recs2[0].flags(), 163);
        assert_eq!(recs2[0].tid(), 0);
        assert_eq!(recs2[0].pos(), 24647210);
    }

    #[test]
    fn test_align_multiple_reads() {
        if !have_refs() {
            return;
        }

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
    fn test_get_record() {
        if !have_refs() {
            return;
        }

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
    fn test_mmap_crasher() {
        if !have_refs() {
            return;
        }

        // this read caused a segfault when there was a bug in the mmap reference loader impl
        let settings = StarSettings::new(DEFAULT_REF_2);
        let reference = StarReference::load(settings).unwrap();
        let mut aligner = reference.get_aligner();

        let read = b"GACAAGCCTGGCCAACATGGTGAAAACCTGTCTCTACTAAAAAAAAAAAAAAATACAAAGATTAGCCGGGTGTGGTGGCAGGCACCTGTAATCCCAGC";
        let qual = b"GGGGAGIGGAAGGGGIIGIGGGIGGGGIIIGGGGGGIIIIIGGGIIIIIIIGGAGAGGAAAGGGAGA.AA<AAGG.GAG.AAA.AAGG<GGGGGA..G";
        let name = b"gatactaga";
        let res = aligner.align_read(name, read, qual);
        assert!(res.len() > 0);
        println!("{:?}", res);
    }

    #[test]
    fn test_write_bam() {
        if !have_refs() {
            return;
        }

        let res = {
            let settings = StarSettings::new(DEFAULT_REF_1);
            let reference = StarReference::load(settings).unwrap();
            let mut aligner = reference.get_aligner();

            let mut out =
                bam::Writer::from_path(&"test/test.bam", &reference.header(), bam::Format::BAM)
                    .unwrap();
            let read = b"GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC";
            let qual = b"GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG";
            let name = b"gatactaga";
            let res = aligner.align_read(name, read, qual);
            for record in res.iter() {
                out.write(&record).unwrap();
            }

            res
            // this scope closes the BAM writing
        };

        let bam_wrapped = bam::Reader::from_path(&"test/test.bam");
        match bam_wrapped {
            Ok(v) => println!("working with version: {:?}", v),
            Err(e) => println!("error parsing header: {:?}", e),
        }
        // The reading portion is commented out because it's failing for unknown reasons, but the BAM
        // file can be read with samtools and pysam.
        let mut bam = bam::Reader::from_path(&"test/test.bam").unwrap();
        let read_records = bam.records();
        let mut i = 0;
        for r in read_records {
            let record = r.unwrap();
            assert_eq!(record.pos(), res[i].pos());
            i += 1;
        }
    }
}
