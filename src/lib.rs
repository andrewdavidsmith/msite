/* MIT License
 *
 * Copyright (c) 2023-2024 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

//! An abstraction for the methylaton level data associated with an
//! individual site in a reference genome, where the data comes from
//! counts of molecules indicating either methylated or unmethylated
//! observations within sequenced reads.
//!
//! Current tasks:
//!
//! - [ ] Replace the `chrom` instance variable with an index that
//!       points to the chrom.
//! - [ ] Allow for a dictionary to be used when the names of
//!       chromosomes are needed.
//! - [ ] Replace the type of `n_reads` with a template so they can be
//!       smaller if desired.
//! - [ ] Replace the `context` with an index to contexts.
//! - [ ] Replace the `meth` variable with an integer value and have
//!       the fractional value calculated when needed.

use std::cmp::max;
use std::cmp::Ordering;
use std::error::Error;
use std::str;

#[derive(Debug, PartialEq, PartialOrd)]
pub struct MSite {
    pub chrom: Vec<u8>,
    pub pos: u64,
    pub strand: char,
    pub context: Vec<u8>,
    pub meth: f64,
    pub n_reads: u64,
}

impl Eq for MSite {}

impl Ord for MSite {
    fn cmp(&self, other: &MSite) -> Ordering {
        self.partial_cmp(&other).unwrap()
    }
}

impl MSite {
    pub fn new() -> MSite {
        MSite {
            chrom: Vec::new(),
            pos: 0,
            strand: '+',
            context: Vec::new(),
            meth: 0.0,
            n_reads: 0,
        }
    }

    pub fn build(s: &String) -> Result<MSite, Box<dyn Error>> {
        let mut parts = s.split_whitespace();

        let chrom = match parts.next() {
            Some(part) => part.to_string().into_bytes(),
            None => return Err("failed to extract chrom".into()),
        };
        let pos = match parts.next() {
            Some(part) => part.parse::<u64>().unwrap(),
            None => return Err("failed to extract pos".into()),
        };
        let strand = match parts.next() {
            Some(part) => part.parse::<char>().unwrap(),
            None => return Err("failed to extract strand".into()),
        };
        let context = match parts.next() {
            Some(part) => part.to_string().into_bytes(),
            None => return Err("failed to extract context".into()),
        };
        let meth = match parts.next() {
            Some(part) => part.parse::<f64>().unwrap(),
            None => return Err("failed to extract meth".into()),
        };
        let n_reads = match parts.next() {
            Some(part) => part.parse::<u64>().unwrap(),
            None => return Err("failed to extract n_reads".into()),
        };

        Ok(MSite {
            chrom,
            pos,
            strand,
            context,
            meth,
            n_reads,
        })
    }
    pub fn n_meth(&self) -> u64 {
        return ((self.n_reads as f64) * self.meth).round() as u64;
    }
    pub fn n_umeth(&self) -> u64 {
        return self.n_reads - self.n_meth();
    }
    pub fn is_cpg(&self) -> bool {
        self.context.len() >= 3
            && self.context[0] == b'C'
            && self.context[1] == b'p'
            && self.context[2] == b'G'
    }
    pub fn is_ccg(&self) -> bool {
        self.context.len() >= 3
            && self.context[0] == b'C'
            && self.context[1] == b'C'
            && self.context[2] == b'G'
    }
    pub fn is_cxg(&self) -> bool {
        self.context.len() >= 3
            && self.context[0] == b'C'
            && self.context[1] == b'X'
            && self.context[2] == b'G'
    }
    pub fn is_chh(&self) -> bool {
        self.context.len() >= 3
            && self.context[0] == b'C'
            && self.context[1] == b'H'
            && self.context[2] == b'H'
    }
    pub fn is_mate_of(&self, second: &MSite) -> bool {
        (self.pos + 1 == second.pos)
            && (self.is_cpg() && second.is_cpg())
            && (self.strand == '+' && second.strand == '-')
    }
    pub fn add(&mut self, other: &MSite) {
        if !self.is_mutated() && other.is_mutated() {
            self.context.push(b'x');
        }
        let total_c_reads = self.n_meth() + other.n_meth();
        self.n_reads += other.n_reads;
        self.meth = (total_c_reads as f64) / (max(1, self.n_reads) as f64);
    }
    pub fn is_mutated(&self) -> bool {
        match self.context.len() {
            0 => false,
            n => self.context[n - 1] == b'x',
        }
    }
    pub fn set_unmutated(&mut self) {
        if self.is_mutated() {
            self.context.pop();
        }
    }
}

impl std::fmt::Display for MSite {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // this method of rounding is to simulate the C++ default
        // formatting of floating point numbers as closely as possible
        const DIGITER: f64 = 1000000.0;
        let m = (self.meth * DIGITER).round() / DIGITER;
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            str::from_utf8(&self.chrom).unwrap(),
            self.pos,
            self.strand,
            str::from_utf8(&self.context).unwrap(),
            m,
            self.n_reads
        )
    }
}
