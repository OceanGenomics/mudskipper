use std::convert::TryFrom;

use rust_htslib::bam::{record::Aux, record::CigarString, HeaderView, Read, Reader, Record};

pub struct BAMQueryRecord {
    is_paired: bool,
    first: Vec<Record>,
    second: Vec<Record>,
}

impl BAMQueryRecord {
    //XXX: Not sure why we can't access the fields directly instead of going through these
    pub fn is_paired(&self) -> bool {
        self.is_paired
    }

    pub fn get_first(&self) -> &Vec<Record> {
        &self.first
    }

    pub fn get_second(&self) -> &Vec<Record> {
        &self.second
    }
}

pub struct BAMQueryRecordReader {
    bam_reader: Reader,
    header: HeaderView,
    last_qname: String,
    record_list: Vec<Record>,
    supp_list: Vec<Record>,
}

impl BAMQueryRecordReader {

    /// Create a new BAMQueryRecordReader that wraps htslib Reader for incremental reading of a BAM / SAM file
    ///
    /// # Arguments
    ///
    /// * `bam_filename` - Path to the BAM file.
    /// * `thread_num` - Optional number of threads to use for reading the BAM file.
    ///
    /// # Panics
    ///
    /// Panics if file isn't in BAM / SAM format or if no non-supplementary records are found in the BAM file.
    ///
    /// # Examples
    ///
    /// ```
    /// use query_bam_records::BAMQueryRecordReader;
    ///
    /// let reader = BAMQueryRecordReader::new("file.bam", Some(4));
    /// ```
    pub fn new(bam_filename: &String, thread_num: Option<usize>) -> BAMQueryRecordReader {
       let mut hts_reader = Reader::from_path(bam_filename).unwrap();
        if let Some(n_threads) = thread_num {
            hts_reader.set_threads(n_threads).expect(&format!("Failed to set number of threads to {}", n_threads))
        }
        let header = hts_reader.header().to_owned();

        // Check to make sure we have a valid BAM / SAM file that contains at least one non-supplementary record
        // by iterating through the records until we find a "representative" record.
        let mut last_qname = String::from("");
        let mut representative_reads: Vec<Record> = Vec::new();
        let mut supplementary_reads: Vec<Record> = Vec::new();
        let mut record = Record::new();

        while let Some(res) = hts_reader.read(&mut record) {
            res.expect("Failed to parse BAM record");
            last_qname = String::from_utf8(record.qname().to_vec()).unwrap();

            // Save records and break when you find the first non-supplementary record
            if record.flags() & 0x800 != 0 {
                supplementary_reads.push(record.to_owned());
            } else {
                representative_reads.push(record.to_owned());
                break;
            }
        }

        if representative_reads.is_empty() {
            panic!("Could not find any non-supplementary records in the BAM file!");
        }

        BAMQueryRecordReader {
            bam_reader: hts_reader,
            header,
            last_qname,
            record_list: representative_reads,
            supp_list: supplementary_reads,
        }
    }

    pub fn get_header(&self) -> &HeaderView {
        &self.header
    }

    /// Create a new Record instance populating it with specific properties of sa_tag (chimeric alignment)
    ///
    /// # Arguments
    ///
    /// * `sa_tag` - SA tag string containing chimeric alignment information.
    ///
    /// # Returns
    ///
    /// A vector containing the parsed `Record` instances representing the chimeric alignments.
    ///
    /// # Examples
    ///
    /// ```
    /// use query_bam_records::BAMQueryRecordReader;
    ///
    /// let reader = BAMQueryRecordReader::new("file.bam", Some(4));
    /// let sa_tag = "1,100,+,10M2D30M,40,60;2,200,-,40M2I20M,35,50;";
    /// let records = reader.get_records_from_sa_tag(sa_tag);
    /// ```

    //XXX: Why is this not used? Can we delete it?
    #[allow(dead_code)]
    fn get_records_from_sa_tag(&self, sa_tag: &str) -> Vec<Record> {
        let mut sa_records: Vec<Record> = Vec::new();
        for aln_str in sa_tag.split(';') {
            if !aln_str.is_empty() {
                let mut brecord = Record::new();
                let tag_vec: Vec<&str> = aln_str.split(',').collect();

                // Parse the fourth value as a byte sequence to create a CigarString instance
                // If parsing fails, panic with an error message
                let cigar: CigarString = CigarString::try_from(tag_vec[3].as_bytes()).expect("Unable to parse cigar string.");

                // Relevant properties of the Record instance (brecord) are set correctly based on the information extracted from the SA tag substrings.
                // If an error occurs it was cause the program to panic
                brecord.set("".as_bytes(), Some(&cigar), "".as_bytes(), "".as_bytes());
                brecord.set_tid(self.header.tid(tag_vec[0].as_bytes()).expect("Cannot find tid for SA alignment!") as i32);
                brecord.set_pos(tag_vec[1].parse::<i64>().expect("Cannot parse position for SA alignment!") - 1);

                // Determines the orientation of the alignment (forward or reverse strand).
                if tag_vec[2].eq("-") {
                    brecord.set_reverse();
                } else {
                    brecord.unset_reverse();
                }

                // Set the mapping quality, If an error occurs it will cause the program to panic.
                brecord.set_mapq(tag_vec[4].parse::<u8>().expect("Cannot parse MAPQ for SA alignment!"));
                sa_records.push(brecord);
            }
        }
        sa_records
    }

    /// Parse the first supplementary alignment of an SA tag into a `Record` instance
    ///
    /// # Arguments
    ///
    /// * `sa_tag` - The SA tag string containing chimeric alignment information.
    ///
    /// # Returns
    ///
    /// A `Record` instance representing the first supplementary alignment, or an empty `Record`, if no alignment information is present
    ///
    /// # Panics
    ///
    /// This function panics if parsing of the SA tag or setting record properties fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use query_bam_record::BAMQueryRecordReader;
    ///
    /// let reader = BAMQueryRecordReader::new("file.bam", Some(4));
    /// let sa_tag = "1,100,+,10M2D30M,40,60;2,200,-,40M2I20M,35,50;";
    /// let primary_record = reader.get_primary_record_of_sa_tag(sa_tag);
    /// ```

    //XXX: According to the BAM spec, SA alignments are in arbitrary order, so I'm not sure why Ehsan calls this the "primary record"
    fn get_primary_record_of_sa_tag(&self, sa_tag: &str) -> Record {
        let mut record = Record::new();
        if let Some(aln_str) = sa_tag.split(';').filter(|s| !s.is_empty()).next() {
            let tag_vec: Vec<&str> = aln_str.split(',').collect();

            // Parse the fourth value as a byte sequence to create a CigarString instance
            // If parsing fails, panic with an error message
            let cigar: CigarString = CigarString::try_from(tag_vec[3].as_bytes()).expect("Unable to parse cigar string.");

            // Relevant properties of the Record instance (record) are set correctly based on the information extracted from the SA tag substrings.
            // If an error occurs it will cause the program to panic.
            record.set("".as_bytes(), Some(&cigar), "".as_bytes(), "".as_bytes());
            record.set_tid(self.header.tid(tag_vec[0].as_bytes()).expect("Cannot find tid for SA alignment!") as i32);
            record.set_pos(tag_vec[1].parse::<i64>().expect("Cannot parse position for SA alignment!") - 1);

            // Determines the orientation of the alignment (forward or reverse strand).
            if tag_vec[2].eq("-") {
                record.set_reverse();
            } else {
                record.unset_reverse();
            }

            // Set the mapping quality, If an error occurs it will cause the program to panic.
            record.set_mapq(tag_vec[4].parse::<u8>().expect("Cannot parse MAPQ for SA alignment!"));
        }
        record
    }

    /// Groups the alignment records into query records based on their primary and supplementary alignments for data organization and optimization.
    ///
    /// # Returns
    ///
    /// A `Result` containing a vector of `BAMQueryRecord` instances if successful, or an error message as a `String` if there is an issue with pairing the alignments.
    ///
    /// # Example
    ///
    /// ```
    /// use query_bam_records::{BAMQueryRecordReader, BAMQueryRecord};
    ///
    /// let reader = BAMQueryRecordReader::new("file.bam", Some(4));
    /// let record_groups = reader.group_records();
    ///
    /// match record_groups {
    ///     Ok(groups) => {
    ///         for group in groups {
    ///             // Process the query records
    ///         }
    ///     }
    ///     Err(error) => {
    ///         println!("Error occurred while grouping records: {}", error);
    ///     }
    /// }
    /// ```

    fn group_records(&self) -> Result<Vec<BAMQueryRecord>, String> {
        let mut record_groups: Vec<BAMQueryRecord> = Vec::new();

        // Find the "primary" record of the SA tag, or panic if there is no SA tag in the record
        let mut primary_of_supp: Vec<Record> = Vec::new();
        for supp in self.supp_list.iter() {
            if let Ok(Aux::String(sa_tag)) = supp.aux(b"SA") {
                primary_of_supp.push(self.get_primary_record_of_sa_tag(sa_tag))
            } else {
                //XXX: Is it really appropriate to panic here, and not return an Error()?
                panic!(
                    "Error reading SA tag for query {}",
                    String::from_utf8(supp.qname().to_vec()).expect("cannot find the qname!")
                );
            }
        }

        // Keep track if a particular supplementary alignment has been assigned to a primary alignment.
        let mut supp_assigned: Vec<bool> = vec![false; self.supp_list.len()];
        let mut i = 0;
        while i < self.record_list.len() {
            let mut first_vec: Vec<Record> = vec![];
            let mut second_vec: Vec<Record> = vec![];
            if !self.record_list[i].is_paired() {
                first_vec.push(self.record_list[i].to_owned());
                // Search for possible matching of a supplementary alignment with the alignment in first_vec.
                for j in 0..supp_assigned.len() {
                    if !supp_assigned[j] && BAMQueryRecordReader::records_equal(&self.record_list[i], &primary_of_supp[j]) {
                        first_vec.push(self.supp_list[j].to_owned());
                        supp_assigned[j] = true;
                    }
                }
                record_groups.push(BAMQueryRecord {
                    is_paired: false,
                    first: first_vec,
                    second: second_vec,
                });
                i += 1;
            } else if self.record_list[i].is_paired() && self.record_list.len() % 2 == 0 {
                first_vec.push(self.record_list[i].to_owned());
                second_vec.push(self.record_list[i + 1].to_owned());
                // search for possible matching of a supplementary alignment with the alignment in first_vec or second_vec
                for j in 0..supp_assigned.len() {
                    if !supp_assigned[j] {
                        if BAMQueryRecordReader::records_equal(&self.record_list[i], &primary_of_supp[j]) {
                            first_vec.push(self.supp_list[j].to_owned());
                            supp_assigned[j] = true;
                        } else if BAMQueryRecordReader::records_equal(&self.record_list[i + 1], &primary_of_supp[j]) {
                            second_vec.push(self.supp_list[j].to_owned());
                            supp_assigned[j] = true;
                        }
                    }
                }
                record_groups.push(BAMQueryRecord {
                    is_paired: true,
                    first: first_vec,
                    second: second_vec,
                });
                i += 2;
            } else {
                return Err(self.last_qname.to_string());
            }
        }
        Ok(record_groups)
    }

    fn records_equal(a: &Record, b: &Record) -> bool {
        a.tid() == b.tid() && a.pos() == b.pos() && a.is_reverse() == b.is_reverse() && a.cigar() == b.cigar()
    }


    /// Groups the alignment records into query records based on their alignment type and matching supplementary alignments, while also handling skipped queries.
    ///
    /// # Returns
    ///
    /// A vector of `BAMQueryRecord` instances representing the grouped query records.
    ///
    /// # Examples
    ///
    /// ```
    /// use my_crate::{BAMQueryRecordReader, BAMQueryRecord};
    ///
    /// let reader = BAMQueryRecordReader::new("file.bam", Some(4));
    /// let record_groups = reader.group_records_skip();
    ///
    /// for group in record_groups {
    ///     // Process the query records
    /// }
    /// ```

    fn group_records_skip(&self) -> Vec<BAMQueryRecord> {
        let mut record_groups: Vec<BAMQueryRecord> = Vec::new();
        // find primary alignment of each supplementary alignment
        let mut primary_of_supp: Vec<Record> = Vec::new();
        for supp in self.supp_list.iter() {
            if let Ok(Aux::String(sa_tag)) = supp.aux(b"SA") {
                primary_of_supp.push(self.get_primary_record_of_sa_tag(sa_tag))
            } else {
                panic!(
                    "Error reading SA tag for query {}",
                    String::from_utf8(supp.qname().to_vec()).expect("cannot find the qname!")
                );
            }
        }
        // Initialize a boolean vector keep track if supplementary alignment has been assigned to a primary alignment
        let mut supp_assigned: Vec<bool> = vec![false; self.supp_list.len()];
        let mut i = 0;
        while i < self.record_list.len() {
            let mut first_vec: Vec<Record> = vec![];
            let mut second_vec: Vec<Record> = vec![];
            if !self.record_list[i].is_paired() {
                first_vec.push(self.record_list[i].to_owned());
                // Search for possible matching of a supplementary alignment with the alignment in first_vec.
                for j in 0..supp_assigned.len() {
                    if !supp_assigned[j] && BAMQueryRecordReader::records_equal(&self.record_list[i], &primary_of_supp[j]) {
                        first_vec.push(self.supp_list[j].to_owned());
                        supp_assigned[j] = true;
                    }
                }
                record_groups.push(BAMQueryRecord {
                    is_paired: false,
                    first: first_vec,
                    second: second_vec,
                });
                i += 1;
            } else if self.record_list[i].is_paired() && self.record_list.len() % 2 == 0 {
                first_vec.push(self.record_list[i].to_owned());
                second_vec.push(self.record_list[i + 1].to_owned());

                // search for possible matching of a supplementary alignment with the alignment in first_vec or second_vec
                for j in 0..supp_assigned.len() {
                    if !supp_assigned[j] {
                        if BAMQueryRecordReader::records_equal(&self.record_list[i], &primary_of_supp[j]) {
                            first_vec.push(self.supp_list[j].to_owned());
                            supp_assigned[j] = true;
                        } else if BAMQueryRecordReader::records_equal(&self.record_list[i + 1], &primary_of_supp[j]) {
                            second_vec.push(self.supp_list[j].to_owned());
                            supp_assigned[j] = true;
                        }
                    }
                }
                record_groups.push(BAMQueryRecord {
                    is_paired: true,
                    first: first_vec,
                    second: second_vec,
                });
                i += 2;
            } else {
                log::warn!("Skipping query: {}", self.last_qname);
                break;
            }
        }
        record_groups
    }

    /// Retrieves the next set of query records from the BAM reader and handle the retrieval and grouping of BAM records based on their query names.
    ///
    /// # Returns
    ///
    /// - `Ok(Some(query_records))`: If there are more query records available, returns the vector of query records.
    /// - `Ok(None)`: If there are no more query records available.
    /// - `Err(error)`: If there is an error while reading the BAM record or grouping the records.
    ///
    /// # Examples
    ///
    /// ```
    /// use my_crate::{BAMQueryRecordReader, BAMQueryRecord};
    ///
    /// let mut reader = BAMQueryRecordReader::new("file.bam", Some(4));
    ///
    /// match reader.get_next_query_records() {
    ///     Ok(Some(query_records)) => {
    ///         for group in query_records {
    ///             // Process the query records
    ///         }
    ///     }
    ///     Ok(None) => {
    ///         println!("No more query records");
    ///     }
    ///     Err(error) => {
    ///         println!("Error occurred: {}", error);
    ///     }
    /// }
    /// ```

    pub fn get_next_query_records(&mut self) -> Result<Option<Vec<BAMQueryRecord>>, String> {
        let mut brecord = Record::new();
        let _query_records: Vec<BAMQueryRecord>;
        while let Some(res) = self.bam_reader.read(&mut brecord) {
            res.expect("Failed to parse BAM record");

            // Extract the query name
            let qname = String::from_utf8(brecord.qname().to_vec()).unwrap();
            if qname != self.last_qname {
                // process the previous group
                match self.group_records() {
                    Ok(query_records) => {
                        self.last_qname = qname;
                        self.record_list.clear();
                        self.supp_list.clear();
                        if brecord.flags() & 0x800 != 0 {
                            self.supp_list.push(brecord.to_owned());
                        } else {
                            self.record_list.push(brecord.to_owned());
                        }
                        return Ok(Some(query_records));
                    }
                    Err(e) => return Err(e),
                };
            } else {
                // The current record belongs to the same query
                // Add the current record to the appropriate list based on the "supplementary" flag
                if brecord.flags() & 0x800 != 0 {
                    self.supp_list.push(brecord.to_owned());
                } else {
                    self.record_list.push(brecord.to_owned());
                }
            }
        }

        match self.group_records() {
            Ok(query_records) => {
                // reset the record and supplementary lists
                self.record_list.clear();
                self.supp_list.clear();
                if query_records.is_empty() {
                    Ok(None)
                } else {
                    Ok(Some(query_records))
                }
            }
            Err(e) => Err(e),
        }
    }

    /// Retrieves the next set of query records from the BAM reader, while skipping queries that cannot be grouped.
    ///
    /// # Returns
    ///
    /// - `Some(query_records)`: If there are more query records available, returns the vector of query records.
    /// - `None`: If there are no more query records available.
    ///
    /// # Examples
    ///
    /// ```
    /// use my_crate::{BAMQueryRecordReader, BAMQueryRecord};
    ///
    /// let mut reader = BAMQueryRecordReader::new("file.bam", Some(4));
    ///
    /// while let Some(query_records) = reader.get_next_query_records_skip() {
    ///     for group in query_records {
    ///         // Process the query records
    ///     }
    /// }
    /// ```

    pub fn get_next_query_records_skip(&mut self) -> Option<Vec<BAMQueryRecord>> {
        let mut brecord = Record::new();
        let query_records: Vec<BAMQueryRecord>;

        // Iterate through the BAM file and read each record
        while let Some(res) = self.bam_reader.read(&mut brecord) {
            res.expect("Failed to parse BAM record");
            let qname = String::from_utf8(brecord.qname().to_vec()).unwrap();
            if qname != self.last_qname {
                // A new query has started, so process the previous group of query records
                query_records = self.group_records_skip();

                // Reset to prepare for processing a new query.
                self.last_qname = qname;
                self.record_list.clear();
                self.supp_list.clear();

                // Add the current record.
                if brecord.flags() & 0x800 != 0 {
                    self.supp_list.push(brecord.to_owned());
                } else {
                    self.record_list.push(brecord.to_owned());
                }
                return Some(query_records);
            } else {
                if brecord.flags() & 0x800 != 0 {
                    self.supp_list.push(brecord.to_owned());
                } else {
                    self.record_list.push(brecord.to_owned());
                }
            }
        }
        // Process the last group of query records and reset members for the next iteration, if it happens
        query_records = self.group_records_skip();
        self.record_list.clear();
        self.supp_list.clear();
        if query_records.is_empty() {
            None
        } else {
            Some(query_records)
        }
    }
}