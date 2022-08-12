use bio_types::strand::Strand;
use mudskipper::convert::{convert_cigar_overhang, TranscriptInfo};
use rust_htslib::bam::record::{CigarString, Cigar};

#[test]
fn test_overhang_cigar() {
    // 100M
    let cigar1 = CigarString(vec![Cigar::Match(100)]);
    // 70M2D28M
    let cigar2 = CigarString(vec![Cigar::Match(70), Cigar::Del(2), Cigar::Match(28)]);
    // 9M1D80M1D9M
    let cigar3 = CigarString(vec![Cigar::Match(9), Cigar::Del(1), Cigar::Match(80), Cigar::Del(1), Cigar::Match(9)]);
    // 7M1D93M
    let cigar4 = CigarString(vec![Cigar::Match(7), Cigar::Del(1), Cigar::Match(93)]);
    // 93M1D7M
    let cigar5 = CigarString(vec![Cigar::Match(93), Cigar::Del(1), Cigar::Match(7)]);
    // 7M1I92M
    let cigar6 = CigarString(vec![Cigar::Match(7), Cigar::Ins(1), Cigar::Match(92)]);
    // 92M1I7M
    let cigar7 = CigarString(vec![Cigar::Match(92), Cigar::Ins(1), Cigar::Match(7)]);
    // 3S97M
    let cigar8 = CigarString(vec![Cigar::SoftClip(3), Cigar::Match(97)]);
    // 97M3S
    let cigar9 = CigarString(vec![Cigar::Match(97), Cigar::SoftClip(3)]);
    
    let tinfo = TranscriptInfo {
        pos: 1,
        strand: Strand::Forward,
        left_overhang: 10,
        right_overhang: 10,
    };

    let converted_cigar1 = convert_cigar_overhang(&cigar1, &tinfo);
    let converted_cigar2 = convert_cigar_overhang(&cigar2, &tinfo);
    let converted_cigar3 = convert_cigar_overhang(&cigar3, &tinfo);
    let converted_cigar4 = convert_cigar_overhang(&cigar4, &tinfo);
    let converted_cigar5 = convert_cigar_overhang(&cigar5, &tinfo);
    let converted_cigar6 = convert_cigar_overhang(&cigar6, &tinfo);
    let converted_cigar7 = convert_cigar_overhang(&cigar7, &tinfo);
    let converted_cigar8 = convert_cigar_overhang(&cigar8, &tinfo);
    let converted_cigar9 = convert_cigar_overhang(&cigar9, &tinfo);

    assert_eq!(format!("{}", converted_cigar1), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar2), "10S60M2D18M10S");
    assert_eq!(format!("{}", converted_cigar3), "10S78M10S");
    assert_eq!(format!("{}", converted_cigar4), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar5), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar6), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar7), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar8), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar9), "10S80M10S");
}