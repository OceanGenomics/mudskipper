use coitrees::COITree;
use rust_htslib::bam::record;
use std::collections::{HashMap, HashSet};
use bio::alphabets::dna::revcomp;

use rust_htslib::bam::HeaderView;
// use rust_htslib::bam::{Format, Header, Read, Reader, Writer, header};
use fnv::FnvHashMap;

extern crate bio_types;
use bio_types::strand::Strand;

use crate::annotation;
use annotation::ExonNode;
use crate::query_bam_records::{BAMQueryRecord};

use log::{error, warn, debug};

/// Given ranges of genomic bases where a spliced alignment is mapped against, returns a map of transcripts that are covered by the projected alignments.
/// For each transcript, the position and orientation of the alignment is also returned.
///
/// # Arguments
///
/// * `tree` - INPUT. The COITree containing exon coordinates of the target chromosome
/// * `ranges` - INPUT. The ranges of genomic bases for the input spliced alignment
pub fn find_tid(tree: &COITree<ExonNode, u32>, ranges: &Vec<(i32, i32)>) -> HashMap<i32, (i32, Strand)> {
    let mut tid_set: HashSet<i32> = HashSet::new(); // stores all the transcripts that are overlapping ALL the ranges of the spliced alignment
    let mut tid_pos: HashMap<i32, (i32, Strand)> = HashMap::new(); // final (tid, pos, strand) to return
    for (i, range) in ranges.iter().enumerate() {
        let range_cov = tree.coverage(range.0, range.1);
        let mut tid_curr: HashSet<i32> = HashSet::new(); // stores all the transcript ids that are covered by the current range 
        debug!("range coverage: {:?}", range_cov);
        if range_cov.0 != 0 || range_cov.1 != 0 {
            tree.query(range.0, range.1, |node| {
                debug!("query for {} {} => {}", range.0, range.1, node.metadata);
                // TODO: for now we are dropping overhanging alignments. This can be improved.
                if node.metadata.start <= range.0 && node.metadata.end >= range.1 {
                    // add to tid_curr
                    if  ranges.len() == 1 || // if there is only a single range, so no need to check splicing boundaries. Otherwise, obey splicing!
                    (i == 0 && range.1 == node.metadata.end) ||
                    (i == ranges.len() - 1 && range.0 == node.metadata.start) || 
                    (i > 0 && i < ranges.len() - 1 && range.0 == node.metadata.start && range.1 == node.metadata.end) {
                        debug!("keeping");
                        tid_curr.insert(node.metadata.tid);
                    }
                    // keep track of the alignment position and strand
                    if i == 0 && node.metadata.strand == Strand::Forward { // first range of the spliced alignment 
                        debug!("inserting => tid:{}, strand:{}, pos:{} = start:{} - tpos_start:{}", 
                            node.metadata.tid, 
                            node.metadata.strand, 
                            node.metadata.start - node.metadata.tpos_start, 
                            node.metadata.start, 
                            node.metadata.tpos_start);
                        tid_pos.insert(node.metadata.tid, (node.metadata.start - node.metadata.tpos_start, Strand::Forward));
                    } else if i == ranges.len() - 1 && node.metadata.strand == Strand::Reverse {
                        debug!("inserting => tid:{}, strand:{}, pos:{} = end:{} - tpos_start:{} + 1", 
                            node.metadata.tid, 
                            node.metadata.strand, 
                            node.metadata.end + node.metadata.tpos_start + 1, 
                            node.metadata.end, 
                            node.metadata.tpos_start);
                        tid_pos.insert(node.metadata.tid, (node.metadata.end + node.metadata.tpos_start + 1, Strand::Reverse));
                    }
                }
            });
            debug!("found coverage: {:?}", range_cov);
        }
        if i ==0 {
            tid_set = tid_curr;
        } else {
            tid_set = &tid_set & &tid_curr; // only transcripts that overlap ALL the ranges will be reported, hence the intersection
        }
    }
    // 
    let mut tid_to_remove: HashSet<i32> = HashSet::new();
    for (k, _v) in &tid_pos {
        if tid_set.contains(k) == false {
            tid_to_remove.insert(*k);
        }
    }
    for tid in tid_to_remove {
        tid_pos.remove(&tid);
    }
    for (k, _v) in &tid_pos {
        debug!("overlapping tid: {}", k);
    }
    return tid_pos;
}

pub fn find_tids_paired(
    tree: &COITree<ExonNode, u32>,
    read_pos1: &i32,
    cigar1: &record::CigarStringView,
    new_cigar1: &mut record::CigarString,
    len1: &mut i32,
    read_pos2: &i32,
    cigar2: &record::CigarStringView,
    new_cigar2: &mut record::CigarString,
    len2: &mut i32,
    long_softclip: &mut bool,
    max_softlen: &usize,
) -> HashMap<i32, ((i32, Strand), (i32, Strand))> {
    let mut tid_pos: HashMap<i32, ((i32, Strand), (i32, Strand))> = HashMap::new();
    debug!("read1: {} {}", read_pos1, cigar1);
    let ranges1 = find_ranges_single(read_pos1, cigar1, new_cigar1, len1, long_softclip, max_softlen);
    let tids1 = find_tid(&tree, &ranges1);
    debug!("read2: {} {}", read_pos2, cigar2);
    let ranges2 = find_ranges_single(read_pos2, cigar2, new_cigar2, len2, long_softclip, max_softlen);
    let tids2 = find_tid(&tree, &ranges2);
    for (tid, pos1) in tids1.iter() {
        if tids2.contains_key(tid) {
            tid_pos.insert(*tid, (*pos1, tids2[tid]));
        }
    }
    return tid_pos;
}

/// Returns a vector of ranges in exons covered by the input alignment. 
/// Each range is a tuple in format of (start, end) and both coordinates are inclusive.
///
/// # Arguments
///
/// * `read_pos` - INPUT. The alignment position of the input BAM record
/// * `cigar` - INPUT. The CIGAR string view of the BAM record
/// * `new_cigar_view` - OUTPUT. The CIGAR string for the alignment projected to the transcriptome
/// * `ref_len` - OUTPUT. The number of reference bases covered by the input alignment
/// * `long_softclip` - OUTPUT. Whether the total number of soft-clipped bases is above the threshold
/// * `max_softlen` - INPUT. The maximum number of soft-clipped bases allowed
pub fn find_ranges_single(
    read_pos: &i32,
    cigar: &record::CigarStringView,
    new_cigar: &mut record::CigarString,
    ref_len: &mut i32,
    long_softclip: &mut bool,
    max_softlen: &usize,
) -> Vec<(i32, i32)> {
    let mut ranges = Vec::new();
    let mut curr_range: (i32, i32) = (-1, -1);

    // This variable declares whether a new exon is reached
    // or the next cigar item corresponds to the current exon positions.
    let mut new_range: bool = true;

    // Set the current position to the base right before the beginning of the alignment
    // debug!("bam_pos: {}", *read_pos);
    let mut curr_pos = *read_pos - 1;
    // debug!("curr_pos: {}", curr_pos);
    let mut new_cigar_vec: Vec<record::Cigar> = Vec::<record::Cigar>::new();
    *ref_len = 0;
    for cigar_item in cigar.iter() {
        debug!("cigar item: {} {}", cigar_item.char(), cigar_item.len());
        let cigar_char = cigar_item.char();
        let cigar_len = cigar_item.len();
        // TODO: Do we need support hard-clip?
        match cigar_char {
            'M' | '=' | 'X' | 'I' | 'D' => {
                if new_range {
                    curr_range = (curr_pos + 1, -1);
                    new_range = false;
                }
                
                if cigar_char != 'I' {
                    // debug!("curr_pos: {}", curr_pos);
                    curr_pos = curr_pos + cigar_len as i32;
                    // debug!("ref_len:{} + cigar_len:{} = {}", *ref_len, cigar_len, *ref_len + cigar_len as i32);
                    *ref_len = *ref_len + cigar_len as i32;
                }

                // debug!("curr_pos: {}", curr_pos);
                curr_range.1 = curr_pos;
                
                // update the new cigar
                if new_cigar_vec.len() == 0 || cigar_char != new_cigar_vec.last().unwrap().char() {
                    new_cigar_vec.push(*cigar_item);
                } else {
                    let new_cigar_len = cigar_len + new_cigar_vec.last().unwrap().len();
                    new_cigar_vec.pop();
                    match cigar_char {
                        'M' => new_cigar_vec.push(record::Cigar::Match(new_cigar_len)),
                        '=' => new_cigar_vec.push(record::Cigar::Equal(new_cigar_len)),
                        'X' => new_cigar_vec.push(record::Cigar::Diff(new_cigar_len)),
                        'I' => new_cigar_vec.push(record::Cigar::Ins(new_cigar_len)),
                        'D' => new_cigar_vec.push(record::Cigar::Del(new_cigar_len)),
                        _ => warn!("Unexpected cigar item: {}", new_cigar_len),
                    }
                }
            }
            // Observing N means that the rest of
            // the cigar belongs to another exon
            'N' => {
                new_range = true;
                ranges.push(curr_range);
                debug!("RANGE {} {}", curr_range.0, curr_range.1);
                curr_pos = curr_pos + cigar_len as i32;
                // debug!("curr_pos: {}", curr_pos);
                // debug!("len:{} + cigar_len:{} = {}", *ref_len, cigar_len, *ref_len + cigar_len as i32);
                *ref_len = *ref_len + cigar_len as i32;
            }
            // Soft clipped bases don't consume reference bases
            'S' => {
                new_cigar_vec.push(*cigar_item);
                if cigar_len > *max_softlen as u32 {
                    *long_softclip = true;
                }
            }
            _ => warn!("Unexpected cigar char! {}", cigar_char),
        }
        // debug!("{} {}", curr_range.0, curr_range.1);
    }
    ranges.push(curr_range);
    debug!("RANGE {} {}", curr_range.0, curr_range.1);
    debug!("REF_LEN {}", *ref_len);
    *new_cigar = record::CigarString(new_cigar_vec);
    return ranges;
}

pub fn convert_paired_end(
    bam_record1: &record::Record,
    bam_record2: &record::Record,
    header_view: &HeaderView,
    transcripts: &Vec<String>,
    txp_lengths: &Vec<i32>,
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    max_softlen: &usize,
) -> Vec<record::Record> {
    let mut converted_records: Vec<record::Record> = Vec::new();
    if bam_record1.is_unmapped() && bam_record2.is_unmapped() {
        converted_records.push(bam_record1.clone());
        converted_records.push(bam_record2.clone());
        return converted_records;
    } else if bam_record1.is_unmapped() && !bam_record2.is_unmapped() {
        let converted_se = convert_single_end(bam_record2, &header_view, transcripts, trees, max_softlen);
        for txp_rec in converted_se.iter() {
            converted_records.push(bam_record1.clone());
            converted_records.push(txp_rec.clone());
        }
        return converted_records;
    } else if !bam_record1.is_unmapped() && bam_record2.is_unmapped() {
        let converted_se = convert_single_end(bam_record1, &header_view, transcripts, trees, max_softlen);
        for txp_rec in converted_se.iter() {
            converted_records.push(txp_rec.clone());
            converted_records.push(bam_record2.clone());
        }
        return converted_records;
    }
    let mut cigar1_new: record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
    let mut cigar2_new: record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
    let mut len1 = 0;
    let mut len2 = 0;
    let mut long_softclip = false;

    let genome_tname = String::from_utf8(header_view.tid2name(bam_record2.tid() as u32).to_vec()).expect("cannot find the tname!");
    if let Some(tree) = trees.get(&genome_tname) {
        let tids = find_tids_paired(
            &tree,
            &(bam_record1.pos() as i32),
            &bam_record1.cigar(),
            &mut cigar1_new,
            &mut len1,
            &(bam_record2.pos() as i32),
            &bam_record2.cigar(),
            &mut cigar2_new,
            &mut len2,
            &mut long_softclip,
            &max_softlen,
        );
        if long_softclip {
            debug!("The softclip length is too long!");
            return converted_records;
        }
        debug!("{}: {}", bam_record1.cigar(), bam_record1.cigar().len());
        debug!("{}: {}", bam_record2.cigar(), bam_record2.cigar().len());
        if tids.len() > 0 {
            for (tid, pos_strand) in tids.iter() {
                let mut first_record_ = bam_record1.clone();
                let mut first_record_cigar = cigar1_new.clone();
                let mut first_record_seq = bam_record1.seq().as_bytes();
                let mut first_record_qual = bam_record1.qual().to_vec();
                let mut second_record_ = bam_record2.clone();
                let mut second_record_cigar = cigar2_new.clone();
                let mut second_record_seq = bam_record2.seq().as_bytes();
                let mut second_record_qual = bam_record2.qual().to_vec();

                let mut first_pos = 0;
                let mut second_pos = 0;
                if pos_strand.0 .1 == Strand::Forward {
                    first_pos = bam_record1.pos() - (pos_strand.0 .0 as i64);
                    debug!("first_pos:{} - pos:{} = {}", bam_record1.pos(), pos_strand.0 .0, first_pos);
                    second_pos = bam_record2.pos() - (pos_strand.1 .0 as i64);
                    debug!("second_pos:{} - pos:{} = {}", bam_record2.pos(), pos_strand.1 .0, second_pos);
                } else if pos_strand.0 .1 == Strand::Reverse {
                    first_pos = (pos_strand.0 .0 as i64) - bam_record1.pos() - len1 as i64;
                    debug!(
                        "pos:{} - first_pos:{} - len1:{} = {}",
                        pos_strand.0 .0,
                        bam_record1.pos(),
                        len1,
                        first_pos
                    );
                    second_pos = (pos_strand.1 .0 as i64) - bam_record2.pos() - len2 as i64;
                    debug!(
                        "pos:{} - second_pos:{} - len2:{} = {}",
                        pos_strand.1 .0,
                        bam_record2.pos(),
                        len2,
                        second_pos
                    );
                    if bam_record1.is_reverse() {
                        debug!("bam_record1 is reversed");
                        first_record_.unset_reverse();
                        first_record_.set_mate_reverse();
                    } else {
                        debug!("bam_record1 is not reversed");
                        first_record_.set_reverse();
                        first_record_.unset_mate_reverse();
                    }
                    if bam_record2.is_reverse() {
                        debug!("second_record is reversed");
                        second_record_.unset_reverse();
                        second_record_.set_mate_reverse();
                    } else {
                        debug!("second_record is not reversed");
                        second_record_.set_reverse();
                        second_record_.unset_mate_reverse();
                    }
                    // reverse the first record
                    let mut cigar_new_rev: Vec<record::Cigar> = Vec::<record::Cigar>::new();
                    for cigar_item in first_record_cigar.into_view(0).iter().rev() {
                        cigar_new_rev.push(*cigar_item);
                    }
                    first_record_cigar = record::CigarString(cigar_new_rev);
                    first_record_seq = revcomp(first_record_seq);
                    first_record_qual.reverse();
                    // reverse the second record
                    let mut cigar_new_rev: Vec<record::Cigar> = Vec::<record::Cigar>::new();
                    for cigar_item in second_record_cigar.into_view(0).iter().rev() {
                        cigar_new_rev.push(*cigar_item);
                    }
                    second_record_cigar = record::CigarString(cigar_new_rev);
                    second_record_seq = revcomp(second_record_seq);
                    second_record_qual.reverse();
                }
                let first_read_len: i64 = bam_record1.seq().len() as i64;
                let second_read_len: i64 = bam_record2.seq().len() as i64;
                let first_length: i64;
                let second_length: i64;
                if pos_strand.0 .1 == Strand::Forward {
                    first_length = if !bam_record1.is_reverse() {
                        second_pos - first_pos + second_read_len
                    } else {
                        second_pos - first_pos + first_read_len
                    };
                    second_length = if !bam_record2.is_reverse() {
                        first_pos - second_pos - first_read_len
                    } else {
                        first_pos - second_pos - second_read_len
                    };
                } else {
                    first_length = if !bam_record1.is_reverse() {
                        second_pos - first_pos - first_read_len
                    } else {
                        second_pos - first_pos - second_read_len
                    };
                    second_length = if !bam_record2.is_reverse() {
                        first_pos - second_pos + second_read_len
                    } else {
                        first_pos - second_pos + first_read_len
                    };
                }
                debug!("tid:{} {}", tid, transcripts[*tid as usize]);
                debug!(
                    "first_pos:{} second_pos:{} len1:{} len2:{} first_length:{} second_length:{}",
                    first_pos, second_pos, first_read_len, second_read_len, first_length, second_length
                );
                debug!("bam_record1.is_reverse():{}", bam_record1.is_reverse());
                debug!("record.is_reverse():{}", bam_record2.is_reverse());
                first_record_.set(bam_record1.qname(), Some(&first_record_cigar), &first_record_seq, &first_record_qual);
                first_record_.set_tid(*tid);
                first_record_.set_mtid(*tid);

                first_record_.set_pos(first_pos);
                first_record_.set_mpos(second_pos);
                first_record_.set_insert_size(first_length);
                // first_record_.remove_aux("AS".as_bytes());

                second_record_.set(bam_record2.qname(), Some(&second_record_cigar), &second_record_seq, &second_record_qual);
                second_record_.set_tid(*tid);
                second_record_.set_mtid(*tid);

                second_record_.set_pos(second_pos);
                second_record_.set_mpos(first_pos);
                second_record_.set_insert_size(second_length);
                // second_record_.remove_aux("AS".as_bytes());

                if first_pos > second_pos {
                    if first_pos + first_read_len > txp_lengths[*tid as usize].into() || second_pos < 0 {
                        continue;
                    }
                    // output_bam.write(&second_record_).unwrap();
                    // output_bam.write(&first_record_).unwrap();
                    converted_records.push(second_record_);
                    converted_records.push(first_record_);
                } else {
                    if second_pos + second_read_len > txp_lengths[*tid as usize].into() || first_pos < 0 {
                        continue;
                    }
                    // output_bam.write(&first_record_).unwrap();
                    // output_bam.write(&second_record_).unwrap();
                    converted_records.push(first_record_);
                    converted_records.push(second_record_);
                }
            }
        } else {
            // missed_read = missed_read + 1;
            debug!("missed read!")
        }
    } else {
        // log for unannotated splicing junction
        debug!("unannotated splicing junction observed!")
    }
    return converted_records;
}

pub fn convert_single_end(
    bam_record: &record::Record,
    header_view: &HeaderView,
    transcripts: &Vec<String>,
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    max_softlen: &usize,
) -> Vec<record::Record> {
    let mut converted_records: Vec<record::Record> = Vec::new();
    if bam_record.is_unmapped() { // no need to convert, just return the unmapped record
        converted_records.push(bam_record.clone());
        return converted_records;
    }

    let mut cigar_new: record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
    let mut ref_len = 0;
    let mut long_softclip = false;

    let ranges = find_ranges_single(
        &(bam_record.pos() as i32),
        &bam_record.cigar(),
        &mut cigar_new,
        &mut ref_len,
        &mut long_softclip,
        &max_softlen,
    );
    let genome_tname = String::from_utf8(header_view.tid2name(bam_record.tid() as u32).to_vec()).expect("cannot find the tname!");
    if let Some(tree) = trees.get(&genome_tname) {
        let tids = find_tid(&tree, &ranges);
        if long_softclip {
            debug!("The softclip length is too long!");
            return converted_records;
        }
        if tids.len() > 0 {
            for (tid, pos_strand) in tids.iter() {
                debug!("tid:{} {}", tid, transcripts[*tid as usize]);

                let mut record_ = bam_record.clone();
                let mut record_cigar = cigar_new.clone();
                let mut record_seq = bam_record.seq().as_bytes();
                let mut record_qual = bam_record.qual().to_vec();
                let mut pos = 0;
                if pos_strand.1 == Strand::Forward {
                    pos = bam_record.pos() - (pos_strand.0 as i64);
                    debug!("bam_pos:{} - pos:{} = {}", bam_record.pos(), pos_strand.0, pos);
                } else if pos_strand.1 == Strand::Reverse {
                    pos = (pos_strand.0 as i64) - bam_record.pos() - ref_len as i64;
                    debug!("pos:{} - bam_pos:{} - len1:{} = {}", pos_strand.0, bam_record.pos(), ref_len, pos);
                    if bam_record.is_reverse() {
                        record_.unset_reverse();
                    } else {
                        record_.set_reverse();
                    }
                    // reverse the record
                    let mut cigar_new_rev: Vec<record::Cigar> = Vec::<record::Cigar>::new();
                    for cigar_item in record_cigar.into_view(0).iter().rev() {
                        cigar_new_rev.push(*cigar_item);
                    }
                    record_cigar = record::CigarString(cigar_new_rev);
                    record_seq = revcomp(record_seq);
                    record_qual.reverse();
                }

                record_.set(bam_record.qname(), Some(&record_cigar), &record_seq, &record_qual);
                record_.set_tid(*tid);
                record_.set_pos(pos);

                // output_bam.write(&record_).unwrap();
                converted_records.push(record_);
            }
        } else {
            // missed_read = missed_read + 1;
            debug!("missed read!")
        }
    } else {
        // log for unannotated splicing junction
        debug!("unannotated splicing junction observed!")
    }
    return converted_records;
}

pub fn convert_query_bam_records(
    qrecord: &BAMQueryRecord,
    header_view: &HeaderView,
    transcripts: &Vec<String>,
    txp_lengths: &Vec<i32>,
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    max_softlen: &usize,
    required_tags: &Vec<&str>,
) -> Vec<record::Record> {
    let mut converted_records: Vec<record::Record> = Vec::new();
    let first_mate = qrecord.get_first();
    let second_mate = qrecord.get_second();

    // drop the record if it contains chimeric alignments (those with supplementary bit set in FLAF, i.e. 0x800)
    // TODO: improve on this by allowing user to choose what to do: https://github.com/OceanGenomics/mudskipper/issues/2
    if first_mate.len() > 1 || second_mate.len() > 1 {
        return converted_records;
    }

    if qrecord.is_paired() {
        debug!("qname1: {}    qname2: {}", String::from_utf8(first_mate[0].qname().to_vec()).unwrap(), String::from_utf8(second_mate[0].qname().to_vec()).unwrap());
        // check for required tags
        for tag in required_tags.iter() {
            if first_mate[0].aux(tag.as_bytes()).is_err() {
                error!("Could not find {} tag for first mate of read {}", tag, String::from_utf8(first_mate[0].qname().to_vec()).unwrap());
                panic!("Some required tags do not exist!");
            }
            if second_mate[0].aux(tag.as_bytes()).is_err() {
                error!("Could not find {} tag for second mate of read {}", tag, String::from_utf8(second_mate[0].qname().to_vec()).unwrap());
                panic!("Some required tags do not exist!");
            }
        }
        let mut txp_records = convert_paired_end(&first_mate[0], &second_mate[0], header_view, transcripts, txp_lengths, trees, max_softlen);
        converted_records.append(&mut txp_records);
    } else {
        debug!("qname: {}", String::from_utf8(first_mate[0].qname().to_vec()).unwrap());
        // check for required tags
        for tag in required_tags.iter() {
            if first_mate[0].aux(tag.as_bytes()).is_err() {
                error!("Could not find {} tag for read {}", tag, String::from_utf8(first_mate[0].qname().to_vec()).unwrap());
                panic!("Some required tags do not exist!");
            }
        }
        let mut txp_records = convert_single_end(&first_mate[0], header_view, transcripts, trees, max_softlen);
        converted_records.append(&mut txp_records);
    }
    return converted_records;
}