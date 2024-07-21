use std::ffi::CStr;
use std::fs::File;
use std::io::BufWriter;
use hmmer_rs::*;
use log::*;
use std::io::Write;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hmmsearch_on_file() {
        let hmms = Hmm::read_hmms_from_path(std::path::Path::new(
            "/Users/jianshuzhao/hmm_bacteria/TIGR03263.HMM",
        ))
        .unwrap();
        assert!(!hmms.is_empty(), "No HMMs were read from the file");
    
        let hmm = &hmms[0]; // Assuming there is at least one HMM
    
        println!("HMM name: {}", hmm.name());
    
        let mut hmmsearch = HmmerPipeline::new(hmm);
        let hmmsearch_result = hmmsearch.run_hmm_on_file(
            hmm,
            std::path::Path::new("/Users/jianshuzhao/Github/hmmsearch_rs/data/test03.faa"),
        );
    
        println!("HMMsearch result: {:?}", hmmsearch_result);
    
        // Safe checking
        unsafe {
            let th = &*hmmsearch_result.c_th;
            if th.nreported > 0 && !th.hit.is_null() {
                for i in 0..th.nreported as isize {
                    let hit = &*th.hit.offset(i);
                    if !hit.dcl.is_null() {
                        let first_hit_name = CStr::from_ptr(hit.name).to_string_lossy();
                        let first_hit_score = hit.score;
                        let first_domain = &*hit.dcl.offset(0);
                        let first_domain_score = first_domain.bitscore;
                        let first_domain_evalue = first_domain.lnP.exp() * (*hmmsearch_result.c_pli).Z;

                        writeln!(output_file, "{}\t{:.4}\t{}\t{}\t{:.4}\t{:.2e}",
                                 hmm.name(),
                                 first_hit_score,
                                 th.nreported,
                                 first_hit_name,
                                 first_domain_score,
                                 first_domain_evalue
                        ).expect("Failed to write hit details");
                    }
                }
            }
        }
    }
}
