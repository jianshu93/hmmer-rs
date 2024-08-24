use log::*;
use std::ffi::CStr;
use std::ffi::CString;

pub struct Hmm {
    pub c_hmm: *mut libhmmer_sys_2::P7_HMM,
}

impl Hmm {
    pub fn read_hmms_from_path(path: &std::path::Path) -> Result<Vec<Hmm>, &'static str> {
        // char          errbuf[eslERRBUFSIZE];
        #[allow(unused_mut)]
        let mut errbuf = CString::new(vec![1; libhmmer_sys_2::eslERRBUFSIZE as usize])
            .unwrap()
            .into_raw();

        let hmmfile = CString::new(path.to_string_lossy().as_bytes())
            .unwrap()
            .into_raw();

        // Open file
        let mut hfp: *mut libhmmer_sys_2::P7_HMMFILE = std::ptr::null_mut();
        let status1 = unsafe {
            libhmmer_sys_2::p7_hmmfile_Open(hmmfile, std::ptr::null_mut(), &mut hfp, errbuf)
        };
        // eslOK = 0
        if status1 != 0 {
            error!("Error in initial reading of HMM file");
            return Err("Error in initial reading of HMM file");
        }
        debug!("HMM file opened successfully");

        // Read each HMM
        let mut hmms = Vec::new();
        loop {
            // ESL_ALPHABET *abc     = NULL;	/* alphabet (set from the HMM file)*/
            // Set to NULL to not force alphabet
            let mut abc: *mut libhmmer_sys_2::ESL_ALPHABET = std::ptr::null_mut();
            // P7_HMM       *hmm     = NULL;
            let mut hmm: *mut libhmmer_sys_2::P7_HMM = std::ptr::null_mut();

            let status2 = unsafe { libhmmer_sys_2::p7_hmmfile_Read(hfp, &mut abc, &mut hmm) };
            if status2 == libhmmer_sys_2::eslEOF as i32 {
                debug!("EOF reached");
                break;
            } else if status2 != 0 {
                error!("Error in reading HMM from opened file");
                return Err("Error in reading HMM from opened file");
            }
            debug!("HMM read successfully");
            hmms.push(Hmm { c_hmm: hmm });
        }

        // retake pointer to free memory
        unsafe {
            let _ = CString::from_raw(errbuf);
            let _ = CString::from_raw(hmmfile);
        }

        Ok(hmms)
    }

    pub fn c_alphabet(&self) -> *const libhmmer_sys_2::ESL_ALPHABET {
        unsafe { (*self.c_hmm).abc }
    }

    pub fn name(&self) -> String {
        return unsafe {
            CStr::from_ptr((*self.c_hmm).name)
                .to_string_lossy()
                .to_string()
        };
    }

    pub fn length(&self) -> u32 {
        unsafe { (*self.c_hmm).M as u32 }
    }

    pub fn acc(&self) -> String {
        let my_acc = unsafe { (*self.c_hmm).acc };
        if my_acc.is_null() {
            "".to_string() // Otherwise we get a segfault
        } else {
            return unsafe { CStr::from_ptr(my_acc).to_string_lossy().to_string() };
        }
    }

    pub fn desc(&self) -> String {
        let my_desc = unsafe { (*self.c_hmm).desc };
        if my_desc.is_null() {
            "".to_string() // Otherwise we get a segfault
        } else {
            return unsafe { CStr::from_ptr(my_desc).to_string_lossy().to_string() };
        }
    }
}

impl Drop for Hmm {
    fn drop(&mut self) {
        unsafe {
            libhmmer_sys_2::p7_hmm_Destroy(self.c_hmm);
        }
    }
}
