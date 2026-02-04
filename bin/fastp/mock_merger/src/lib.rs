use rayon::prelude::*;
use pyo3::prelude::*;
use std::fs::{OpenOptions, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use log::{info, warn};
use std::time::Instant;
use std::sync::OnceLock;

#[derive(Debug)]
struct FastqRead {
    header: String,
    sequence: String,
    quality: String,
}

static LOGGER_INIT: OnceLock<()> = OnceLock::new();

fn init_logger() {
    LOGGER_INIT.get_or_init(|| {
        env_logger::builder().format_timestamp(None).init();
    });
}

fn read_chunk(reader: &mut dyn BufRead, count: usize) -> std::io::Result<Vec<FastqRead>> {
    let mut records = Vec::with_capacity(count);
    let mut header = String::new();
    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();

    for _ in 0..count {
        header.clear();
        if reader.read_line(&mut header)? == 0 { break; }

        sequence.clear();
        if reader.read_line(&mut sequence)? == 0 { break; }

        plus.clear();
        if reader.read_line(&mut plus)? == 0 { break; }

        quality.clear();
        if reader.read_line(&mut quality)? == 0 { break; }

        records.push(FastqRead {
            header: header.trim_end().to_string(),
            sequence: sequence.trim_end().to_string(),
            quality: quality.trim_end().to_string(),
        });
    }

    Ok(records)
}

fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(|c| match c {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => 'A'
    }).collect()
}

fn merge_reads(r1: &FastqRead, r2: &FastqRead, dist: usize) -> Vec<u8> {
    format!(
        "{} mock_merged_{}_{}\n{}{}{}\n+\n{}{}{}\n",
        r1.header, r1.sequence.len(), r2.sequence.len(),
        r1.sequence, "N".repeat(dist), reverse_complement(&r2.sequence),
        r1.quality, "#".repeat(dist), r2.quality,
    ).into_bytes()
}

#[pyfunction]
fn mock_merge_by_chunks(
    in_fq1: &str,
    in_fq2: &str,
    dist: usize,
    chunk_size: usize,
    out_path: &str,
) -> PyResult<()> {
    init_logger();

    let start = Instant::now();

    info!("Opening FASTQ files: '{}' and '{}'", in_fq1, in_fq2);
    let mut reader1 = BufReader::new(File::open(in_fq1)?);
    let mut reader2 = BufReader::new(File::open(in_fq2)?);
    let mut writer = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(out_path)?
    );
    info!("Output will be written to '{}'", out_path);

    let mut chunk_idx = 0;
    loop {
        info!("Reading chunk {}", chunk_idx);
        let read1_chunk = read_chunk(&mut reader1, chunk_size)?;
        let read2_chunk = read_chunk(&mut reader2, chunk_size)?;

        if read1_chunk.is_empty() || read2_chunk.is_empty() {
            info!("Reached EOF (empty chunk). Stopping.");
            break;
        }

        if read1_chunk.len() != read2_chunk.len() {
            warn!("Mismatched reads in chunk {}: {} vs {}", chunk_idx, read1_chunk.len(), read2_chunk.len());
        }

        info!("Merging {} read pairs", read1_chunk.len().min(read2_chunk.len()));
        let merged: Vec<Vec<u8>> = read1_chunk.par_iter()
            .zip(read2_chunk.par_iter())
            .map(|(a, b)| merge_reads(a, b, dist))
            .collect();

        for record in merged {
            writer.write_all(&record)?;
        }

        chunk_idx += 1;
    }

    let duration = start.elapsed();
    info!("Merge complete in {:.2?}. Output written to '{}'", duration, out_path);
    writer.flush()?;
    Ok(())
}

#[pymodule]
fn mock_merge_rs(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(mock_merge_by_chunks, m)?)?;
    Ok(())
}
