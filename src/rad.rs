use std::io::{stdout, BufReader, BufWriter, Cursor, Seek, SeekFrom, Write};

struct RadRecord {
    pub fw: bool,
    pub qname: String,
    pub tid_list: Vec<i32>
}

// A number of lines in this function is borrowed from bam2rad funciton at
// https://github.com/COMBINE-lab/alevin-fry/blob/master/libradicl/src/convert.rs
pub fn read_bamfile(records: &Vec<RadRecord>) {
    let mut data = Cursor::new(vec![]);

    for record in records.iter() {
        tid_list = record.tid_list;
        if !tid_list.is_empty() {
            assert!(!tid_list.is_empty(), "Trying to write empty tid_list");
            let na = tid_list.len();
            data.write_all(&(na as u32).to_le_bytes()).unwrap();
            // write bc
            // data.write_all(&(bc as u32).to_le_bytes()).unwrap();
            // write umi
            // data.write_all(&(umi as u32).to_le_bytes()).unwrap();
            // write tid list
            for t in tid_list.iter() {
                data.write_all(&t.to_le_bytes()).unwrap();
            }
        }

    }
}