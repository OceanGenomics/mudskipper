use std::convert::TryFrom;

use rust_htslib::bam::{record::Aux, record::CigarString, HeaderView, Read, Reader, Record, ext::BamRecordExtensions};

pub struct BAMQueryRecord {
    is_paired: bool,
    first: Vec<Record>,
    second: Vec<Record>,
}

impl BAMQueryRecord {
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
    /// Returns a reference to the header of the BAM file.
    ///
    /// # Examples
    ///
    /// ```
    /// use my_crate::BAMQueryRecordReader;
    ///
    /// let reader = BAMQueryRecordReader::new("file.bam", Some(4));
    /// let header = reader.get_header();
    /// ```

    pub fn get_header(&self) -> &HeaderView {
        &self.header
    }


    /// Retrieves a list of chimeric alignments from the SA tag.
    ///
    /// # Arguments
    ///
    /// * `sa_tag` - SA tag string containing chimeric alignment information.
    ///
    /// # Returns
    ///
    /// A vector containing the parsed `Record` objects representing the chimeric alignments.
    ///
    /// # Examples
    ///
    /// ```
    /// use my_crate::BAMQueryRecordReader;
    ///
    /// let reader = BAMQueryRecordReader::new("file.bam", Some(4));
    /// let sa_tag = "chr1,100,+,60M;chr2,200,-,40M";
    /// let records = reader.get_records_from_sa_tag(sa_tag);
    /// ```

    #[allow(dead_code)]
    fn get_records_from_sa_tag(&self, sa_tag: &str) -> Vec<Record> {

        // Create an empty vector to store the records
        let mut sa_records: Vec<Record> = Vec::new();

        // Split the SA tag string into substrings using the ';' delimiter
        for aln_str in sa_tag.split(';') {

            // Check if the substring is not empty
            if !aln_str.is_empty() {

                // Create a new Record object
                let mut brecord = Record::new();

                // Split the substring into values using the ',' delimiter and store them in a vector
                let tag_vec: Vec<&str> = aln_str.split(',').collect();


                // Parse the fourth value as a byte sequence to create a CigarString object
                // If parsing fails, panic with an error message
                let cigar: CigarString = CigarString::try_from(tag_vec[3].as_bytes()).expect("Unable to parse cigar string.");

                // Set the properties of the brecord object based on the values in tag_vec
                brecord.set("".as_bytes(), Some(&cigar), "".as_bytes(), "".as_bytes());
                brecord.set_tid(self.header.tid(tag_vec[0].as_bytes()).expect("Cannot find tid for SA alignment!") as i32);
                brecord.set_pos(tag_vec[1].parse::<i64>().expect("Cannot parse position for SA alignment!") - 1);

                // Check if the third value is "-", indicating a reverse flag
                // Set the reverse flag accordingly
                if tag_vec[2].eq("-") {
                    brecord.set_reverse();
                } else {
                    brecord.unset_reverse();
                }

                // Parse the fifth value as an unsigned 8-bit integer to set the mapping quality
                brecord.set_mapq(tag_vec[4].parse::<u8>().expect("Cannot parse MAPQ for SA alignment!"));

                // Add the brecord object to the sa_records vector
                sa_records.push(brecord);
            }
        }

        // Return the vector containing the records
        sa_records
    }

    /// Retrieves the primary alignment record from the SA tag.
    ///
    /// The SA tag is a string containing chimeric alignment information, where the first alignment
    /// represents the primary alignment. This function parses the SA tag and returns the primary
    /// alignment as a `Record` object.
    ///
    /// # Arguments
    ///
    /// * `sa_tag` - The SA tag string containing chimeric alignment information.
    ///
    /// # Returns
    ///
    /// The primary alignment `Record` object.
    ///
    /// # Panics
    ///
    /// This function panics if parsing of the SA tag or setting record properties fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use my_crate::BAMQueryRecordReader;
    ///
    /// let reader = BAMQueryRecordReader::new("file.bam", Some(4));
    /// let sa_tag = "chr1,100,+,60M;chr2,200,-,40M";
    /// let primary_record = reader.get_primary_record_of_sa_tag(sa_tag);
    /// ```

    fn get_primary_record_of_sa_tag(&self, sa_tag: &str) -> Record {

        // Create a new Record object
        let mut brecord = Record::new();

        // Split the SA tag string into substrings using the ';' delimiter
        for aln_str in sa_tag.split(';') {

            // Check if the substring is not empty
            if !aln_str.is_empty() {

                // Split the substring into values using the ',' delimiter and store them in a vector
                let tag_vec: Vec<&str> = aln_str.split(',').collect();

                // Parse the fourth value as a byte sequence to create a CigarString object
                // If parsing fails, panic with an error message
                let cigar: CigarString = CigarString::try_from(tag_vec[3].as_bytes()).expect("Unable to parse cigar string.");

                // Set the properties of the brecord object based on the values in tag_vec
                brecord.set("".as_bytes(), Some(&cigar), "".as_bytes(), "".as_bytes());
                brecord.set_tid(self.header.tid(tag_vec[0].as_bytes()).expect("Cannot find tid for SA alignment!") as i32);
                brecord.set_pos(tag_vec[1].parse::<i64>().expect("Cannot parse position for SA alignment!") - 1);

                // Check if the third value is "-", indicating a reverse flag
                // Set the reverse flag accordingly
                if tag_vec[2].eq("-") {
                    brecord.set_reverse();
                } else {
                    brecord.unset_reverse();
                }

                // Parse the fifth value as an unsigned 8-bit integer to set the mapping quality
                brecord.set_mapq(tag_vec[4].parse::<u8>().expect("Cannot parse MAPQ for SA alignment!"));

                // Return the brecord object as the primary record
                return brecord;
            }
        }
        // If no valid records are found, return the empty brecord object
        brecord
    }

    /// Groups the alignment records into query records based on their primary and supplementary alignments.
    ///
    /// This function takes a list of alignment records and groups them into query records, considering
    /// both single-end and paired-end alignments. Each query record consists of one or two alignment
    /// records representing the primary alignment(s) and any corresponding supplementary alignments.
    /// The grouping is based on equality of the target ID (tid), position (pos), reverse flag, and
    /// CIGAR string between the alignments.
    ///
    /// # Returns
    ///
    /// A `Result` containing a vector of `BAMQueryRecord` objects if successful, or an error message as a `String` if there is an issue with pairing the alignments.
    ///
    /// # Examples
    ///
    /// ```
    /// use my_crate::{BAMQueryRecordReader, BAMQueryRecord};
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

    fn records_equal(a: &Record, b: &Record) -> bool {
        a.tid() == b.tid() && a.pos() == b.pos() && a.is_reverse() == b.is_reverse() && a.cigar() == b.cigar()
    }

    fn group_records(&self) -> Result<Vec<BAMQueryRecord>, String> {
        let mut record_groups: Vec<BAMQueryRecord> = Vec::new();

        // find primary alignment of each supplementary alignment
        let mut primary_of_supp: Vec<Record> = Vec::new();
        for supp in self.supp_list.iter() {
            if let Ok(Aux::String(sa_tag)) = supp.aux(b"SA") {

                // If the SA tag exists and is of type `String`, retrieve the primary record
                primary_of_supp.push(self.get_primary_record_of_sa_tag(sa_tag))
            } else {

                // If the SA tag does not exist or is not of type `String`, panic with an error message
                panic!(
                    "Error reading SA tag for query {}",
                    String::from_utf8(supp.qname().to_vec()).expect("cannot find the qname!")
                );
            }
        }
        // initialize a boolean vector supp_assigned with a length equal to the number of supplementary alignments
        // Each element is set to false initially to keep track of whether a supplementary alignment has been assigned to a primary alignment.
        let mut supp_assigned: Vec<bool> = vec![false; self.supp_list.len()];
        let mut i = 0;
        while i < self.record_list.len() {
            let mut first_vec: Vec<Record> = vec![];
            let mut second_vec: Vec<Record> = vec![];
            if !self.record_list[i].is_paired() {

                // single-end
                // add the primary alignment first
                first_vec.push(self.record_list[i].to_owned());

                // search for possible matching of a supplementary alignment with the alignment in first_vec
                for j in 0..supp_assigned.len() {
                    if !supp_assigned[j] && BAMQueryRecordReader::records_equal(&self.record_list[i], &primary_of_supp[j]) {

                        // println!("ADDING to first_vec");
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
                // paired-end
                // add the primary alignments first
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
                // return the name of the query which missing its mate
                return Err(self.last_qname.to_string());
            }
        }
        Ok(record_groups)
    }

    /// Groups the alignment records into query records based on their alignment type and matching supplementary alignments, while also handling skipped queries.
    ///
    /// This function takes a list of alignment records and groups them into query records, considering both single-end and paired-end alignments. Each query record consists of one or two alignment records representing the primary alignment(s) and any corresponding supplementary alignments. The grouping is based on equality of the target ID (tid), position (pos), reverse flag, and CIGAR string between the alignments.
    ///
    /// # Returns
    ///
    /// A vector of `BAMQueryRecord` objects representing the grouped query records.
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

                // If the SA tag does not exist or is not of type `String`, panic with an error message
                panic!(
                    "Error reading SA tag for query {}",
                    String::from_utf8(supp.qname().to_vec()).expect("cannot find the qname!")
                );
            }
        }
        // is set to false initially to keep track of whether a supplementary alignment has been assigned to a primary alignment
        let mut supp_assigned: Vec<bool> = vec![false; self.supp_list.len()];
        let mut i = 0;
        while i < self.record_list.len() {
            let mut first_vec: Vec<Record> = vec![];
            let mut second_vec: Vec<Record> = vec![];
            if !self.record_list[i].is_paired() {
                // single-end
                // add the primary alignment first
                first_vec.push(self.record_list[i].to_owned());
                // search for possible matching of a supplementary alignment with the alignment in first_vec
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
                // paired-end
                // add the primary alignments first
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

                // print out warning message for skipped query
                log::warn!("Skipping query: {}", self.last_qname);
                break;
            }
        }
        record_groups
    }

    /// Retrieves the next set of query records from the BAM reader.
    ///
    /// This function reads the next BAM record from the reader and groups the records into query records. It returns a `Result` containing either the next set of query records as a vector of `BAMQueryRecord` objects or an error message as a `String`.
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
        // Read the next BAM record using the BAM reader
        // log::debug!("in func");
        while let Some(res) = self.bam_reader.read(&mut brecord) {
            res.expect("Failed to parse BAM record");

            // Extract the query name from the current record
            // log::debug!("in loop!");
            let qname = String::from_utf8(brecord.qname().to_vec()).unwrap();
            if qname != self.last_qname {

                // a new query
                // process the previous group
                match self.group_records() {
                    Ok(query_records) => {

                        // Update the last query name
                        self.last_qname = qname;

                        // Clear the record and supplementary lists
                        self.record_list.clear();
                        self.supp_list.clear();

                        // add the current record
                        if brecord.flags() & 0x800 != 0 {
                            self.supp_list.push(brecord.to_owned());
                        } else {
                            self.record_list.push(brecord.to_owned());
                        }

                        // Return the query records
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
        // process the last record

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
    /// This function reads the next BAM record from the reader and groups the records into query records, skipping queries that cannot be grouped due to missing mate or other conditions. It returns an `Option` containing either the next set of query records as a vector of `BAMQueryRecord` objects or `None` if there are no more query records available.
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

        // log::debug!("in func");
        while let Some(res) = self.bam_reader.read(&mut brecord) {

            // log::debug!("in loop!");
            res.expect("Failed to parse BAM record");
            let qname = String::from_utf8(brecord.qname().to_vec()).unwrap();
            if qname != self.last_qname {

                // a new query
                // process the previous group
                query_records = self.group_records_skip();

                // reset
                self.last_qname = qname;
                self.record_list.clear();
                self.supp_list.clear();

                // add the current record
                if brecord.flags() & 0x800 != 0 {
                    self.supp_list.push(brecord.to_owned());
                } else {
                    self.record_list.push(brecord.to_owned());
                }
                return Some(query_records);
            } else {
                // same query
                if brecord.flags() & 0x800 != 0 {
                    self.supp_list.push(brecord.to_owned());
                } else {
                    self.record_list.push(brecord.to_owned());
                }
            }
        }

        // process the last record
        query_records = self.group_records_skip();

        // reset
        self.record_list.clear();
        self.supp_list.clear();

        if query_records.is_empty() {
            None
        } else {
            Some(query_records)
        }
    }
}
