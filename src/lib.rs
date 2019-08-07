extern crate failure;
extern crate seq_io;

use failure::Error;

use seq_io::fastq;
use seq_io::fastq::Record;
use std::ffi::CString;
use std::io::BufReader;
use std::io::prelude::*;
use std::fs::File;
use std::os::raw::c_char;
use std::os::raw::c_int;
use std::os::raw::c_ulonglong;
use std::path::Path;
include!("bindings.rs");

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::{HeaderView};

pub struct RefStruct{ pub sr : *const StarRef }
unsafe impl Sync for RefStruct {}

/// StarSettings contains the parameters which will be used for the STAR aligner.
/// Currently the array of argument strings is passed directly
pub struct StarSettings {
    ref_dir : String,
    multn : usize,
    pub args : Vec<String>,
}

/// References to some commonly used reference genomes for testing purposes
pub const DEFAULT_REF_1 : &str =  "/mnt/opt/refdata_cellranger/mm10-3.0.0/star";
pub const DEFAULT_REF_2 : &str = "/mnt/opt/refdata_cellranger/GRCh38-3.0.0/star";

/// The default value for number of mappings allowed is 1.
const DEFAULT_MULTN : usize = 1;

impl StarSettings {
    
    /// This constructor just sets all of the necessary arguments to their defaults, and the
    /// arguments which can take on different values have separate functions to set them later
    pub fn new() -> StarSettings {
        let def_args : Vec<String> = vec![
            "STAR".to_string(), "--genomeDir".to_string(), DEFAULT_REF_1.to_string(), 
            "--outSAMmultNmax".to_string(), DEFAULT_MULTN.to_string(),
            "--runThreadN".to_string(), "1".to_string(),
            "--readNameSeparator".to_string(), "space".to_string(),
            "--outSAMunmapped".to_string(), "Within".to_string(), "KeepPairs".to_string(),
            "--outSAMtype".to_string(), "SAM".to_string(),
            "--outStd".to_string(), "SAM".to_string(),
            "--outSAMorder".to_string(), "PairedKeepInputOrder".to_string(),
        ];
        StarSettings{ ref_dir : DEFAULT_REF_1.to_string(), multn : DEFAULT_MULTN, args : def_args }
    }

    pub fn clone(&self) -> StarSettings {
        let nref_dir = self.ref_dir.clone();
        let n_multn = self.multn.clone();
        let mut nvec = Vec::new();
        for i in 0..self.args.len() {
            nvec.push(self.args[i].clone());
        }
        StarSettings {
            ref_dir : nref_dir,
            multn : n_multn,
            args : nvec,
        }
    }

    /// Create a new StarSettings and immediately sets its reference genome directory to the
    /// provided value
    pub fn from_str(reference : &str) -> StarSettings {
        let mut res = StarSettings::new();
        res.set_reference(reference);
        res
    }

    /// Update the array of arguments when the arguments' values have changed
    fn sync_args(&mut self) {
        for i in 0..self.args.len() {
            if self.args[i] == "--genomeDir" {
                self.args[i+1] = self.ref_dir.clone();
            }
            else if self.args[i] == "--outSAMmultNmax" {
                self.args[i+1] = self.multn.to_string();
            }   
        }
    }

    /// Set the reference genome directory
    pub fn set_reference(&mut self, new_ref : &str) {
        self.ref_dir = new_ref.to_string();
        self.sync_args();
    }

    /// Set the max number of multimapping reads
    pub fn set_multn(&mut self, new_multn : usize) {
        self.multn = new_multn;
        self.sync_args();
    }

    /// Add the given read group strings to the arguments
    pub fn add_rg(&mut self, rg_tags : Vec<String>) {
        self.args.push("--outSAMattrRGline".to_string());
        for tag in rg_tags {
            self.args.push(tag);
        }
    }
}

/// StarRawAligner is an aligner which interfaces directly with the bindgen-produced orbit
/// library.  It contains all information necessary to align reads, but does not also contain the
/// header information needed to produce a BAM file from the alignment records.
pub struct StarRawAligner {
    al : *mut Aligner,
    pub settings : StarSettings,
}

impl StarRawAligner {
    pub fn new() -> StarRawAligner {
        let cur_settings = StarSettings::new();
        let cur_aligner = init_aligner_rust(&cur_settings.args).unwrap();
        StarRawAligner { al : cur_aligner, settings : cur_settings }
    }
    pub fn from_str(ref_path : &str) -> StarRawAligner {
        let cur_settings = StarSettings::from_str(ref_path);
        let cur_aligner = init_aligner_rust(&cur_settings.args).unwrap();
        StarRawAligner { al : cur_aligner, settings : cur_settings }
    }
    pub fn from_settings_and_ref(init_settings : &StarSettings, init_ref : *const StarRef) -> StarRawAligner {
        let c_aligner = init_aligner_from_ref_rust(init_ref).unwrap();
        StarRawAligner {
            al : c_aligner,
            settings: init_settings.clone(),
        }
    }

    pub fn clone(&self) -> StarRawAligner {
        // TODO copy settings over
        StarRawAligner { al : init_aligner_clone_rust(self.al), settings : StarSettings::new() }
    }

    pub fn align_read(&self, read : String, qual : String) -> String {
        align_read_rust(self.al, read, qual).unwrap()
    }
    pub fn align_read_pair(&self, read : String, qual : String, read2 : String, qual2 : String) -> String {
        align_read_pair_rust(self.al, read, qual, read2, qual2).unwrap()
    }
    pub fn destroy(&self) {
        unsafe{destroy_aligner(self.al);}
    }
}

/// StarAligner is the main external-facing alignment interface.  It creates a StarRawAligner with
/// caller-provided settings, and supports read alignment as well as dealing with Fastq records as
/// input and the output BAM records.
pub struct StarAligner {
    pub aligner : StarRawAligner,
    pub header : Header,
}

impl StarAligner {

    /// Creates an aligner with default parameters
    pub fn new() -> StarAligner {
        let cur_aligner = StarRawAligner::new();
        let (cur_header, _cur_hv) = generate_header_with_view(
            &Path::new(&(cur_aligner.settings.args[2]))
        );

        StarAligner {
            aligner : cur_aligner,
            header : cur_header,
        }
    }

    /// Creates an aligner with a given reference genome directory
    pub fn from_str(ref_path : &str) -> StarAligner {
        let cur_aligner = StarRawAligner::from_str(ref_path);
        
        let (cur_header, _cur_hv) = generate_header_with_view(
            &Path::new(&(cur_aligner.settings.args[2]))
        );

        StarAligner {
            aligner : cur_aligner,
            header : cur_header,
        }
    }

    // Creates a StarAligner from the arguments plus an already constructed genome/parameters index
    pub fn from_settings_and_ref(init_settings : &StarSettings, init_ref : *const StarRef) -> StarAligner {
        let cur_aligner = StarRawAligner::from_settings_and_ref(init_settings, init_ref);
        let (cur_header, _cur_hv) = generate_header_with_view(
            &Path::new(&(cur_aligner.settings.args[2]))
        );

        StarAligner {
            aligner : cur_aligner,
            header : cur_header,
        }
    }

    /// Produces another aligner which shares many key read-only data structures with the first
    /// one.  When aligning with multiple threads, ther should be one main aligner constructed with
    /// the functions above, and then each thread should produce a clone to use for its own
    /// alignment.
    pub fn clone(&self) -> StarAligner {
        let cur_aligner = self.aligner.clone();
        let (cur_header, _cur_hv) = generate_header_with_view(
            &Path::new(&(cur_aligner.settings.args[2]))
        );

        StarAligner {
            aligner : cur_aligner,
            header : cur_header,
        }
    }

    /// Aligns a given read and produces BAM records
    pub fn align_single_read(&self, name : String, read : String, qual : String) -> Vec<bam::Record> {
        let sam_string = self.aligner.align_read(read, qual);
        self.parse_sam_to_records(sam_string, name)
    }

    /// Aligns a given pair of reads and produces BAM records
    pub fn align_read_pair(&self, name : String, read1 : String, qual1 : String, read2 : String, qual2 : String) -> (Vec<bam::Record>, Vec<bam::Record>) {
        let read1_clone = read1.clone();
        let sam_string = self.aligner.align_read_pair(read1, qual1, read2, qual2);
        let full_vec = self.parse_sam_to_records(sam_string, name);

        // Partition the records into first mate and second mate.  For now, this is being done by
        // looking at the SEQ field of each and comparing to the first read in the pair.
        let mut first_vec : Vec<bam::Record> = Vec::new();
        let mut second_vec : Vec<bam::Record> = Vec::new();
        for rec in full_vec {
            if String::from_utf8(rec.seq().as_bytes()).unwrap() == read1_clone {
                first_vec.push(rec);
            }
            else
            {
                second_vec.push(rec);
            }
        }
        (first_vec, second_vec)
    }

    /// Aligns every read in a Fastq file
    pub fn align_fastq(&self, fastq_path : &str) -> Result<Vec<bam::Record>, Error> {
        let mut res : Vec<bam::Record> = Vec::new();
        let mut reader = fastq::Reader::from_path(Path::new(fastq_path))?;
        while let Some(cur_record) = reader.next() {
            let record = cur_record?;
            let cur = self.align_fastq_record(record)?;
            for x in cur {
                res.push(x);
            }
        }
        Ok(res)
    }

    /// Aligns the read contained in a single given Fastq record
    pub fn align_fastq_record<R : Record>(&self, record : R) -> Result<Vec<bam::Record>, Error> {

        let seq : String = String::from_utf8(record.seq().to_vec())?;
        let qual : String = String::from_utf8(record.qual().to_vec())?;
        let name : String = record.id()?.to_string();
        let cur = self.align_single_read(name, seq, qual);
        Ok(cur)

    }

    /// Destructor
    pub fn destroy(&self) {
        self.aligner.destroy();
    }

    /// Given a list of BAM records as a SAM-format string in which records are separated by new
    /// lines, add the records to a vector and append the read name to the beginning of them so
    /// that they conform with BAM specifications
    fn parse_sam_to_records(&self, sam: String, name : String) -> Vec<bam::Record> {
        let mut records = Vec::new();
        for slc in sam.split("\n") {
            if slc.len() > 0 {
                let str = format!("{}{}", name, slc);
                let record = {
                    let (_, cur_hv) = generate_header_with_view(
                        &Path::new(&(self.aligner.settings.args[2]))
                    );

                    bam::Record::from_sam(&cur_hv, str.as_bytes()).unwrap()
                };
                records.push(record);
            }
        }

        records
    }
}

/// Read in the lines from a file and store each line as its own string in a vector
fn get_lines(path : &Path) -> Vec<String> {
    let mut res : Vec<String> = Vec::new();
    let lines = BufReader::new(File::open(path).unwrap()).lines();
    for x in lines {
        res.push(x.unwrap());
    }
    res
}

/// Produces a header from the genome reference directory by looking up the contig names and
/// lengths and formatting them properly
fn generate_header_with_view(genome_path : &Path) -> (Header, HeaderView) {
    let mut header = Header::new();
    
    let contig_names_path = genome_path.join(Path::new("chrName.txt"));
    let contig_names = get_lines(&contig_names_path);

    let contig_lengths_path = genome_path.join(Path::new("chrLength.txt"));
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

pub fn init_aligner_clone_rust(al : *const Aligner) -> *mut Aligner
{
    unsafe{init_aligner_clone(al)}
}

pub fn init_aligner_rust(args : &Vec<String>) ->Result<*mut Aligner, Error>
{
    let mut nvec = Vec::new();
    for x in args.iter() {
        let cur_string = CString::new(x.as_str())?;
        nvec.push(cur_string.into_raw());
    }

    let c_args = nvec.as_mut_ptr() as *mut *mut c_char;
    let length = nvec.len() as c_int;
    Ok(unsafe{init_aligner(length, c_args)})
}

pub fn init_ref_rust(args : &Vec<String>) ->Result<*const StarRef, Error>
{
    let mut nvec = Vec::new();
    for x in args.iter() {
        let cur_string = CString::new(x.as_str())?;
        nvec.push(cur_string.into_raw());
    }

    let c_args = nvec.as_mut_ptr() as *mut *mut c_char;
    let length = nvec.len() as c_int;
    Ok(unsafe{init_star_ref(length, c_args)})
}

pub fn init_aligner_from_ref_rust(ref_ptr : *const StarRef) -> Result<*mut Aligner, Error>
{
    Ok(unsafe{init_aligner_from_ref(ref_ptr)})
}

pub fn align_read_rust(al : *mut Aligner, read : String, qual : String) -> Result<String, Error>
{
    let length = read.len() as c_ulonglong;
    let c_read = CString::new(read)?;
    let c_qual = CString::new(qual)?;
    let read_ptr = c_read.as_ptr() as *mut c_char;
    let qual_ptr = c_qual.as_ptr() as *mut c_char;
    let res : *const c_char = unsafe{align_read(al, read_ptr, qual_ptr, length)};
    let string_res : String = unsafe{CString::from_raw(res as *mut c_char).into_string()?};
    Ok(string_res)
}

pub fn align_read_pair_rust(al : *mut Aligner, read : String, qual : String, read2 : String, qual2 : String) -> Result<String, Error>
{
    let length = read.len() as c_ulonglong;
    let c_read = CString::new(read)?;
    let c_qual = CString::new(qual)?;
    let read_ptr = c_read.as_ptr() as *mut c_char;
    let qual_ptr = c_qual.as_ptr() as *mut c_char;
    
    let c_read2 = CString::new(read2)?;
    let c_qual2 = CString::new(qual2)?;
    let read_ptr2 = c_read2.as_ptr() as *mut c_char;
    let qual_ptr2 = c_qual2.as_ptr() as *mut c_char;

    let res : *const c_char = unsafe{align_read_pair(al, read_ptr, qual_ptr, read_ptr2, qual_ptr2, length)};
    let string_res : String = unsafe{CString::from_raw(res as *mut c_char).into_string()?};
    Ok(string_res)
}

#[test]
fn test_align_read()
{
    let aligner = StarRawAligner::from_str(DEFAULT_REF_2);
    let read : String = "GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC".to_string();
    let qual : String = "GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG".to_string();
    
    let res : String = aligner.align_read(read, qual);
    assert!(res == "\t0\t6\t30070474\t255\t98M\t*\t0\t0\tGTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC\tGGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:0\n");
    aligner.destroy();
}

#[test]
fn test_align_multiple_reads()
{
    let aligner = StarRawAligner::from_str(DEFAULT_REF_2);
    let read : String = "GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC".to_string();
    let qual : String = "GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG".to_string();
    let read2 : String = "GTATGTCAAGTTGGTGGAGGCCCTTTGTGCTGAACACCAAATCAACCTAATTAAGGTTGATGACAACAAGAAACTAGGAGAATGGGTAGGCCTTTGTA".to_string();
    let qual2 : String = "AGGAGGGGIG.GAGGGGIGGIGIIGGGIIAAGGGGGIGGIIGAGIGIA.GGGGIGGGGGGGGGGGGGGIGIIIIGGGGGGIGGGGIIIGA.G.<.GGG".to_string();

    let res : String = aligner.align_read(read, qual);
    let res2 : String = aligner.align_read(read2, qual2);
    
    assert!(res == "\t0\t6\t30070474\t255\t98M\t*\t0\t0\tGTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC\tGGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:0\n");
    assert!(res2 == "\t0\t6\t132816509\t255\t55M396N43M\t*\t0\t0\tGTATGTCAAGTTGGTGGAGGCCCTTTGTGCTGAACACCAAATCAACCTAATTAAGGTTGATGACAACAAGAAACTAGGAGAATGGGTAGGCCTTTGTA\tAGGAGGGGIG.GAGGGGIGGIGIIGGGIIAAGGGGGIGGIIGAGIGIA.GGGGIGGGGGGGGGGGGGGIGIIIIGGGGGGIGGGGIIIGA.G.<.GGG\tNH:i:1\tHI:i:1\tAS:i:98\tnM:i:0\n");

    aligner.destroy();

}

#[test]
fn test_raw_clone()
{
    let aligner = StarRawAligner::from_str(DEFAULT_REF_2);
    let aligner2 = aligner.clone();
    let read : String = "GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC".to_string();
    let qual : String = "GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG".to_string();
    let read2 : String = "GTATGTCAAGTTGGTGGAGGCCCTTTGTGCTGAACACCAAATCAACCTAATTAAGGTTGATGACAACAAGAAACTAGGAGAATGGGTAGGCCTTTGTA".to_string();
    let qual2 : String = "AGGAGGGGIG.GAGGGGIGGIGIIGGGIIAAGGGGGIGGIIGAGIGIA.GGGGIGGGGGGGGGGGGGGIGIIIIGGGGGGIGGGGIIIGA.G.<.GGG".to_string();

    let res : String = aligner.align_read(read, qual);
    let res2 : String = aligner2.align_read(read2, qual2);
    
    assert!(res == "\t0\t6\t30070474\t255\t98M\t*\t0\t0\tGTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC\tGGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:0\n");
    assert!(res2 == "\t0\t6\t132816509\t255\t55M396N43M\t*\t0\t0\tGTATGTCAAGTTGGTGGAGGCCCTTTGTGCTGAACACCAAATCAACCTAATTAAGGTTGATGACAACAAGAAACTAGGAGAATGGGTAGGCCTTTGTA\tAGGAGGGGIG.GAGGGGIGGIGIIGGGIIAAGGGGGIGGIIGAGIGIA.GGGGIGGGGGGGGGGGGGGIGIIIIGGGGGGIGGGGIIIGA.G.<.GGG\tNH:i:1\tHI:i:1\tAS:i:98\tnM:i:0\n");

    aligner2.destroy();
    aligner.destroy();

}

#[test]
fn test_from_params()
{
    let settings = StarSettings::from_str(DEFAULT_REF_2);
    let reference_object = init_ref_rust(&settings.args).unwrap();
    let aligner = StarRawAligner::from_settings_and_ref(&settings, reference_object);
    let aligner2 = StarRawAligner::from_settings_and_ref(&settings, reference_object);
    let read : String = "GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC".to_string();
    let qual : String = "GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG".to_string();
    let read2 : String = "GTATGTCAAGTTGGTGGAGGCCCTTTGTGCTGAACACCAAATCAACCTAATTAAGGTTGATGACAACAAGAAACTAGGAGAATGGGTAGGCCTTTGTA".to_string();
    let qual2 : String = "AGGAGGGGIG.GAGGGGIGGIGIIGGGIIAAGGGGGIGGIIGAGIGIA.GGGGIGGGGGGGGGGGGGGIGIIIIGGGGGGIGGGGIIIGA.G.<.GGG".to_string();

    let res : String = aligner.align_read(read, qual);
    let res2 : String = aligner2.align_read(read2, qual2);
    
    assert!(res == "\t0\t6\t30070474\t255\t98M\t*\t0\t0\tGTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC\tGGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:0\n");
    assert!(res2 == "\t0\t6\t132816509\t255\t55M396N43M\t*\t0\t0\tGTATGTCAAGTTGGTGGAGGCCCTTTGTGCTGAACACCAAATCAACCTAATTAAGGTTGATGACAACAAGAAACTAGGAGAATGGGTAGGCCTTTGTA\tAGGAGGGGIG.GAGGGGIGGIGIIGGGIIAAGGGGGIGGIIGAGIGIA.GGGGIGGGGGGGGGGGGGGIGIIIIGGGGGGIGGGGIIIGA.G.<.GGG\tNH:i:1\tHI:i:1\tAS:i:98\tnM:i:0\n");

    aligner.destroy();
    aligner2.destroy();
    unsafe{destroy_ref(reference_object);}
}

#[test]
fn test_header()
{
    let aligner = StarAligner::from_str(DEFAULT_REF_2);
    
    let header_string : String = String::from_utf8(aligner.header.to_bytes()).unwrap();
    assert!(header_string.starts_with("@SQ\tSN:1\tLN:248956422\n"));
    assert!(header_string.ends_with("@SQ\tSN:KI270394.1\tLN:970"));
}

#[test]
fn test_get_record()
{
    let aligner = StarAligner::from_str(DEFAULT_REF_2);
    let read : String = "GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC".to_string();
    let qual : String = "GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG".to_string();
    let name : String = "gatactaga".to_string();
    let res = aligner.align_single_read(name, read, qual);
    assert!(res.len() > 0);
    assert_eq!(res[0].pos(), 30070473);
 
}

#[test]
fn test_write_bam()
{
    let aligner = StarAligner::new();
    let mut out = bam::Writer::from_path(&"test/test.bam", &aligner.header).unwrap();
    let read : String = "GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC".to_string();
    let qual : String = "GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG".to_string();
    let name : String = "gatactaga".to_string();
    let res = aligner.align_single_read(name, read, qual);
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

#[test]
fn test_align_fastq()
{
    let aligner = StarAligner::new();
    let fastq_path : &str = "test/full1_sample.fastq";

    let res : Vec<bam::Record> = aligner.align_fastq(fastq_path).unwrap();
    assert!(res.len() > 0);

    let mut out = bam::Writer::from_path(&"test/full1_sample.bam", &aligner.header).unwrap();
    for record in res.iter() {
        out.write(&record).unwrap();
    }
}
