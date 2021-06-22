use bio::io::gff;

pub fn read(ann_file_adr: String) {
    let ann_file_adr_split: Vec<&str> = ann_file_adr.split(".").collect();
    let file_type: &str = ann_file_adr_split.last().copied().unwrap_or("default string");
    println!("{}", file_type);
    let ann_type: gff::GffType = if file_type == "gtf" { gff::GffType::GTF2 } 
                                 else { if file_type == "gff3" { gff::GffType::GFF3 }
                                        else { gff::GffType::GFF2 } };

    let reader = gff::Reader::from_file(ann_file_adr, ann_type);
    let mut n: i32 = 0;
    for record in reader.expect("Error reading file.").records() {
        let rec = record.ok().expect("Error reading record.");
        n += 1;
        print!("\r{} {:?}", n, rec.seqname());
    }
    println!("\nn: {}", n);
}