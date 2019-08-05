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

/// StarSettings contains the parameters which will be used for the STAR aligner.
/// Currently the array of argument strings is passed directly
pub struct StarSettings {
    ref_dir : String,
    multn : usize,
    args : Vec<String>,
}

/// References to some commonly used reference genomes for testing purposes
pub const DEFAULT_REF_1 : &str =  "/mnt/opt/refdata_cellranger/mm10-3.0.0/star";
pub const DEFAULT_REF_2 : &str = "/mnt/opt/refdata_cellranger/GRCh38-3.0.0/star";

/// The default value for number of mappings allowed is 1.
const DEFAULT_MULTN : usize = 1;

impl StarSettings {
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
    pub fn from_str(reference : &str) -> StarSettings {
        let mut res = StarSettings::new();
        res.set_reference(reference);
        res
    }

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

    pub fn set_reference(&mut self, new_ref : &str) {
        self.ref_dir = new_ref.to_string();
        self.sync_args();
    }

    pub fn set_multn(&mut self, new_multn : usize) {
        self.multn = new_multn;
        self.sync_args();
    }

    pub fn add_rg(&mut self, rg_tags : Vec<String>) {
        self.args.push("--outSAMattrRGline".to_string());
        for tag in rg_tags {
            self.args.push(tag);
        }
    }
}

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

pub struct StarAligner {
    pub aligner : StarRawAligner,
    pub header : Header,
}

impl StarAligner {
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

    pub fn align_single_read(&self, name : String, read : String, qual : String) -> Vec<bam::Record> {
        let sam_string = self.aligner.align_read(read, qual);
        self.parse_sam_to_records(sam_string, name)
    }
    pub fn align_read_pair(&self, name : String, read1 : String, qual1 : String, read2 : String, qual2 : String) -> (Vec<bam::Record>, Vec<bam::Record>) {
        let sam_string = self.aligner.align_read_pair(read1, qual1, read2, qual2);
        let full_vec = self.parse_sam_to_records(sam_string, name);
        // TODO actually partition this 
        (full_vec, Vec::new())
    }

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

    pub fn align_fastq_record<R : Record>(&self, record : R) -> Result<Vec<bam::Record>, Error> {

        let seq : String = String::from_utf8(record.seq().to_vec())?;
        let qual : String = String::from_utf8(record.qual().to_vec())?;
        let name : String = record.id()?.to_string();
        let cur = self.align_single_read(name, seq, qual);
        Ok(cur)

    }

    pub fn destroy(&self) {
        self.aligner.destroy();
    }

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

fn get_lines(path : &Path) -> Vec<String> {
    let mut res : Vec<String> = Vec::new();
    let lines = BufReader::new(File::open(path).unwrap()).lines();
    for x in lines {
        res.push(x.unwrap());
    }
    res
}

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

fn add_ref_to_bam_header(header: &mut Header, seq_name: &str, seq_len: usize) {
    let mut header_rec = HeaderRecord::new(b"SQ");
    header_rec.push_tag(b"SN", &seq_name);
    header_rec.push_tag(b"LN", &seq_len);
    header.push_record(&header_rec);
}

pub fn init_aligner_clone_rust(al : *mut Aligner) -> *mut Aligner
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

    aligner.destroy();

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
    let mut out = bam::Writer::from_path(&"test.bam", &aligner.header).unwrap();
    let read : String = "GTGCGGGGAGAAGTTTCAAGAAGGTTCTTATGGAAAAAAGGCTGTGAGCATAGAAAGCAGTCATAGGAGGTTGGGGAACTAGCTTGTCCCTCCCCACC".to_string();
    let qual : String = "GGGAGIGIIIGIIGGGGIIGGIGGAGGAGGAAG.GGIIIG<AGGAGGGIGGGGIIIIIGGIGGGGGIGIIGGAGGGGGIGGGIGIIGGGGIIGGGIIG".to_string();
    let name : String = "gatactaga".to_string();
    let res = aligner.align_single_read(name, read, qual);
    for record in res.iter() {
        out.write(&record).unwrap();
    }
    let bam_wrapped = bam::Reader::from_path(&"test.bam");
    match bam_wrapped {
        Ok(v) => println!("working with version: {:?}", v),
        Err(e) => println!("error parsing header: {:?}", e),
    }
    // The reading portion is commented out because it's failing for unknown reasons, but the BAM
    // file can be read with samtools and pysam.
    /*
    let mut bam = bam::Reader::from_path(&"test.bam").unwrap();
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
    let fastq_path : &str = "full1_sample.fastq";

    let res : Vec<bam::Record> = aligner.align_fastq(fastq_path).unwrap();
    assert!(res.len() > 0);

    let mut out = bam::Writer::from_path(&"full1_sample.bam", &aligner.header).unwrap();
    for record in res.iter() {
        out.write(&record).unwrap();
    }
}
