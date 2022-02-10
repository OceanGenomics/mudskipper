use std::path::Path;
use std::cmp;
use std::hash::Hasher;
use std::collections::hash_map::DefaultHasher;
use std::fs;
use std::str;
use rust_htslib::bam;
use rust_htslib::bam::{HeaderView, record::Aux, Read};
use log::{trace, info};

pub fn depositionify_bam(input_path: &str, output_path: &str, max_mem: u64, nthreads: usize) {
    let mut bam = threaded_bam_reader(input_path, nthreads);
    let header = bam::Header::from_template(bam.header());
    let bam_bytes = fs::metadata(&input_path).unwrap().len();
    let buckets = ((bam_bytes / max_mem) + 1) as u32;

    let output = if input_path == output_path {
        let mut hasher = DefaultHasher::new();
        hasher.write(output_path.as_bytes());
        format!("{}.{}.tmp", output_path, hasher.finish())
    } else {
        output_path.to_string()
    };
    
    info! {"Processing file {} of size {} into {} buckets", input_path, bam_bytes, buckets };
    clean_buckets(); 
    bucketify_bam(&mut bam, &header, buckets);
    process_buckets(&header, &output, buckets, nthreads);
    clean_buckets(); 

    if output != output_path {
        fs::remove_file(output_path).unwrap();
        fs::rename(output, output_path).unwrap();
    }
}

fn bucketify_bam(bam: &mut bam::Reader, header: &bam::Header, nbuckets: u32) {
    fs::create_dir_all(Path::new("buckets")).unwrap();
    let mut writers: Vec<_> = bucket_names("buckets", nbuckets).iter()
        .map(|p| bam::Writer::from_path(p, &header, bam::Format::Bam).unwrap())
        .collect();
    
    for r in bam.records() {
        let record = match r {
            Ok(r) => {
                trace!("Processing record: {:?}", r);
                r
            }
            Err(e) => {
                info! { "failed to read BAM record {:?}", e };
                continue
            }
        };
        //Grab the read name and figure out which bucket it goes into
        let qname = record.qname();
        let mut hasher = DefaultHasher::new();
        hasher.write(qname);
        let hash = hasher.finish();
        let bucket = (hash % (nbuckets as u64)) as usize;
        //XXX: Not sure if this is buffered I/O or not, if so, we'll have to manually buffer
        writers[bucket].write(&record).expect("Failed to write record");
    }
}

fn process_buckets(header: &bam::Header, output_path: &str, nbuckets: u32, nthreads: usize) {
    info! { "Gathering buckets into reordered BAM file at {:?}", output_path }
    let output_path = Path::new(output_path);
    fs::create_dir_all(output_path.parent().unwrap()).expect("failed to create output directory");
    if output_path.exists() {
        fs::remove_file(output_path).expect("output file exists and cannot be removed")
    }
    let mut writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam).unwrap();

    let mut readers: Vec<_> = bucket_names("buckets", nbuckets).iter()
        .map(|b| (b.to_string(), threaded_bam_reader(b, nthreads)))
        .collect();
    
    for (name, reader) in &mut readers {
        info! { "Reading bucket {:?}", name}
        process_bucket(header, reader, &mut writer)
    }
}

fn process_bucket(header: &bam::Header, reader: &mut bam::Reader, writer: &mut bam::Writer) {
    //XXX: avoid re-hashing these records every time
    let mut records: Vec<_> = reader.records().filter_map(|r| r.ok()).collect();
    let header_view = HeaderView::from_header(header);
    records.sort_by_cached_key(|r| sort_key(&header_view, r));
    for record in records {
        writer.write(&record).unwrap();
    }
}

fn bucket_names(prefix: &str, nbuckets: u32) -> Vec<String> {
    (0..nbuckets).map(|i| format!("{}/bucket_{}.bam", prefix, i)).collect()
}

fn sort_key(header: &bam::HeaderView, r: &bam::Record) -> (u64, u64, u64, u64) {
    //We hash both the read name and the position, combine them, and return that as the sort key
    //This to eliminate problems with unusual read names (that could have some meaningful order),
    //and to ensure that multiple alignments from one read don't come out in any particular order
    let mut qhasher = DefaultHasher::new();
    qhasher.write(r.qname());

    //For paired reads, hash out the lowest reference ID and position
    let mut mhasher = DefaultHasher::new();
    mhasher.write_i32(cmp::min(r.tid(), r.mtid()));
    mhasher.write_i32(cmp::max(r.tid(), r.mtid()));
    mhasher.write_i64(cmp::min(r.pos(), r.mpos()));
    mhasher.write_i64(cmp::max(r.pos(), r.mpos()));

    //For chimeric reads, again hash out the lowest reference ID and position
    let mut chash: u64 = 0;
    if let Ok(Aux::String(sa_tag)) = r.aux(b"SA") {
        let mut _hash: u64 = 0;
        for aln_str in sa_tag.split(";") {
            if aln_str.is_empty() == false {
                //info! { "parsing aln_str {:?}", aln_str}
                let tag_vec: Vec<&str> = aln_str.split(",").collect();
                let tid = header.tid(tag_vec[0].as_bytes()).unwrap();
                let pos = tag_vec[1].parse::<i64>().unwrap();
                let strand = tag_vec[2];
                let mut chasher = DefaultHasher::new();
                chasher.write_i32(cmp::min(tid as i32, r.tid()));
                chasher.write_i32(cmp::max(tid as i32, r.tid()));
                chasher.write_i64(cmp::min(pos, r.pos()));
                chasher.write_i64(cmp::max(pos, r.pos()));
                chasher.write(strand.as_bytes());
                chash = chasher.finish();
            }
        }
    };

    //Finally, find the position of the particular read and randomize that
    let mut poshasher = DefaultHasher::new();
    poshasher.write_i64(r.pos());
    (qhasher.finish(), mhasher.finish(), chash, poshasher.finish())
}

fn threaded_bam_reader(path: &str, _nthreads: usize) -> bam::Reader {
    //XXX: There's some sort of terrible threading issue that's causing truncated reads
    bam::Reader::from_path(path).unwrap()

    // let mut reader = bam::Reader::from_path(path).unwrap();
    //     So, disable threading for now
    // if nthreads > 1 {
    //     reader.set_threads((nthreads as usize) - 1).unwrap();
    // } else {
    //     reader.set_threads(1).unwrap();
    // }
    // reader
}

fn clean_buckets() {
	let _ = fs::remove_dir_all(Path::new("buckets")); //XXX: do something safer here
}