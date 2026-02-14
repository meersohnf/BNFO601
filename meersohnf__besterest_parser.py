"""
YOURNAME_besterest_parser.py

A robust BLAST XML parser implementing a finite state machine architecture.
This parser handles BLAST XML output format (outfmt 5) and is version-independent.

Author: YOURNAME
Course: BNFO 601
Date: February 2026

Based on betterer_parser.py but adapted for XML format for improved stability.
"""

from dataclasses import dataclass
from decimal import Decimal
from typing import Tuple, List


# Data classes for structured storage of BLAST results
# These use modern Python dataclass decorator for clean, immutable data structures

@dataclass(frozen=True)
class Hsp:
    """
    Represents a High-Scoring Segment Pair - one local alignment between query and subject.
    HSPs are the fundamental unit of BLAST alignments.
    """
    bit_score: str
    score: str
    evalue: str
    query_from: int
    query_to: int
    hit_from: int
    hit_to: int
    identity: int
    positive: int
    gaps: int
    align_len: int
    qseq: str
    hseq: str
    midline: str

    @property
    def evalue_decimal(self):
        """Convert e-value to Decimal for arbitrary precision"""
        return Decimal(self.evalue)

    @property
    def percent_identity(self):
        """Calculate percent identity of alignment"""
        return (self.identity / self.align_len * 100) if self.align_len > 0 else 0.0


@dataclass(frozen=True)
class Hit:
    """
    Represents one subject sequence that matched the query.
    A hit can contain multiple HSPs if there are multiple local alignment regions.
    """
    hit_num: int
    hit_id: str
    hit_def: str
    hit_accession: str
    hit_len: int
    hsps: Tuple[Hsp, ...]

    @property
    def best_evalue(self):
        """Return the best (lowest) e-value among all HSPs for this hit"""
        if not self.hsps:
            return None
        return min(hsp.evalue_decimal for hsp in self.hsps)


@dataclass(frozen=True)
class QueryResult:
    """
    Represents one query (iteration) and all its hits.
    This is the main output unit - one per query protein.
    """
    iter_num: int
    query_id: str
    query_def: str
    query_len: int
    hits: Tuple[Hit, ...]
    message: str = ""  # Will be "No hits found" if no hits

    @property
    def has_hits(self):
        """Check if query has any hits"""
        return len(self.hits) > 0

    @property
    def hit_count(self):
        """Return number of hits"""
        return len(self.hits)


class BlastXMLParser:
    """
    A finite state machine parser for BLAST XML output.

    This parser is robust and version-independent because XML provides
    a stable, self-describing format unlike plain text BLAST output.

    The parser operates as an FSM with these states:
    - scan_for_iteration: Looking for next query
    - in_iteration: Extracting query metadata
    - scan_for_hit: Looking for next hit or checking for no hits
    - in_hit: Extracting hit metadata
    - scan_for_hsp: Looking for next HSP
    - in_hsp: Extracting HSP alignment data
    """

    def __init__(self):
        """Initialize parser state"""
        self._reset_all_state()
        self.state = self._scan_for_iteration

    def parse(self, filename: str):
        """
        Parse BLAST XML file and yield QueryResult objects.

        This is a generator that yields one QueryResult per query (iteration).
        Usage: for result in parser.parse('blast_output.xml'): ...
        """
        with open(filename, 'r') as ifh:
            for line in ifh:
                line = line.strip()
                result = self.state(line)
                if result is not None:
                    yield result

        # Emit final query if file ends mid-iteration
        final = self._finalize_query()
        if final is not None:
            yield final

    # ==================== STATE METHODS ====================
    # Each method corresponds to one state in the FSM
    # Methods consume one line and optionally emit a QueryResult

    def _scan_for_iteration(self, line: str):
        """State 1: Scanning for start of new query iteration"""
        if line.startswith('<Iteration>'):
            # Starting a new query - finalize any previous one
            result = self._finalize_query()
            self._reset_query_state()
            self.state = self._in_iteration
            return result
        return None

    def _in_iteration(self, line: str):
        """State 2: Inside iteration, extracting query metadata"""
        if line.startswith('<Iteration_iter-num>'):
            self.iter_num = int(self._extract_content(line))
        elif line.startswith('<Iteration_query-ID>'):
            self.query_id = self._extract_content(line)
        elif line.startswith('<Iteration_query-def>'):
            self.query_def = self._extract_content(line)
        elif line.startswith('<Iteration_query-len>'):
            self.query_len = int(self._extract_content(line))
        elif line.startswith('<Iteration_hits>'):
            self.state = self._scan_for_hit
        elif line.startswith('<Iteration_message>'):
            # Capture "No hits found" message
            self.message = self._extract_content(line)
        return None

    def _scan_for_hit(self, line: str):
        """State 3: Scanning for next hit or end of hits"""
        if line.startswith('<Hit>'):
            self._reset_hit_state()
            self.state = self._in_hit
        elif line.startswith('</Iteration_hits>') or line.startswith('</Iteration>'):
            # No more hits for this query
            result = self._finalize_query()
            self.state = self._scan_for_iteration
            return result
        elif 'No hits found' in line or self.message == 'No hits found':
            # Explicit no hits message
            result = self._finalize_query()
            self.state = self._scan_for_iteration
            return result
        return None

    def _in_hit(self, line: str):
        """State 4: Inside hit, extracting hit metadata"""
        if line.startswith('<Hit_num>'):
            self.hit_num = int(self._extract_content(line))
        elif line.startswith('<Hit_id>'):
            self.hit_id = self._extract_content(line)
        elif line.startswith('<Hit_def>'):
            self.hit_def = self._extract_content(line)
        elif line.startswith('<Hit_accession>'):
            self.hit_accession = self._extract_content(line)
        elif line.startswith('<Hit_len>'):
            self.hit_len = int(self._extract_content(line))
        elif line.startswith('<Hit_hsps>'):
            self.state = self._scan_for_hsp
        return None

    def _scan_for_hsp(self, line: str):
        """State 5: Scanning for next HSP or end of HSPs"""
        if line.startswith('<Hsp>'):
            self._reset_hsp_state()
            self.state = self._in_hsp
        elif line.startswith('</Hit_hsps>'):
            # Done with this hit - save it and look for next hit
            self._finalize_hit()
            self.state = self._scan_for_hit
        return None

    def _in_hsp(self, line: str):
        """State 6: Inside HSP, extracting alignment data"""
        if line.startswith('<Hsp_bit-score>'):
            self.hsp_bit_score = self._extract_content(line)
        elif line.startswith('<Hsp_score>'):
            self.hsp_score = self._extract_content(line)
        elif line.startswith('<Hsp_evalue>'):
            self.hsp_evalue = self._extract_content(line)
        elif line.startswith('<Hsp_query-from>'):
            self.hsp_query_from = int(self._extract_content(line))
        elif line.startswith('<Hsp_query-to>'):
            self.hsp_query_to = int(self._extract_content(line))
        elif line.startswith('<Hsp_hit-from>'):
            self.hsp_hit_from = int(self._extract_content(line))
        elif line.startswith('<Hsp_hit-to>'):
            self.hsp_hit_to = int(self._extract_content(line))
        elif line.startswith('<Hsp_identity>'):
            self.hsp_identity = int(self._extract_content(line))
        elif line.startswith('<Hsp_positive>'):
            self.hsp_positive = int(self._extract_content(line))
        elif line.startswith('<Hsp_gaps>'):
            self.hsp_gaps = int(self._extract_content(line))
        elif line.startswith('<Hsp_align-len>'):
            self.hsp_align_len = int(self._extract_content(line))
        elif line.startswith('<Hsp_qseq>'):
            self.hsp_qseq = self._extract_content(line)
        elif line.startswith('<Hsp_hseq>'):
            self.hsp_hseq = self._extract_content(line)
        elif line.startswith('<Hsp_midline>'):
            self.hsp_midline = self._extract_content(line)
        elif line.startswith('</Hsp>'):
            # Done with this HSP - save it and look for next HSP
            self._finalize_hsp()
            self.state = self._scan_for_hsp
        return None

    # ==================== HELPER METHODS ====================

    def _extract_content(self, line: str) -> str:
        """
        Extract content between XML tags.
        Example: "<Tag>content</Tag>" returns "content"
        """
        # Find content between first '>' and last '<'
        start = line.find('>') + 1
        end = line.rfind('<')
        if start > 0 and end > start:
            return line[start:end]
        return ""

    def _finalize_hsp(self):
        """Create HSP object and add to current hit's HSP list"""
        hsp = Hsp(
            bit_score=self.hsp_bit_score,
            score=self.hsp_score,
            evalue=self.hsp_evalue,
            query_from=self.hsp_query_from,
            query_to=self.hsp_query_to,
            hit_from=self.hsp_hit_from,
            hit_to=self.hsp_hit_to,
            identity=self.hsp_identity,
            positive=self.hsp_positive,
            gaps=self.hsp_gaps,
            align_len=self.hsp_align_len,
            qseq=self.hsp_qseq,
            hseq=self.hsp_hseq,
            midline=self.hsp_midline
        )
        self.current_hsps.append(hsp)

    def _finalize_hit(self):
        """Create Hit object and add to current query's hit list"""
        hit = Hit(
            hit_num=self.hit_num,
            hit_id=self.hit_id,
            hit_def=self.hit_def,
            hit_accession=self.hit_accession,
            hit_len=self.hit_len,
            hsps=tuple(self.current_hsps)
        )
        self.current_hits.append(hit)

    def _finalize_query(self):
        """Create QueryResult object for completed query"""
        if not self.query_id:
            return None

        result = QueryResult(
            iter_num=self.iter_num,
            query_id=self.query_id,
            query_def=self.query_def,
            query_len=self.query_len,
            hits=tuple(self.current_hits),
            message=self.message
        )
        return result

    def _reset_all_state(self):
        """Reset all parser state variables"""
        self._reset_query_state()

    def _reset_query_state(self):
        """Reset state for new query"""
        self.iter_num = 0
        self.query_id = ""
        self.query_def = ""
        self.query_len = 0
        self.message = ""
        self.current_hits = []
        self._reset_hit_state()

    def _reset_hit_state(self):
        """Reset state for new hit"""
        self.hit_num = 0
        self.hit_id = ""
        self.hit_def = ""
        self.hit_accession = ""
        self.hit_len = 0
        self.current_hsps = []
        self._reset_hsp_state()

    def _reset_hsp_state(self):
        """Reset state for new HSP"""
        self.hsp_bit_score = ""
        self.hsp_score = ""
        self.hsp_evalue = ""
        self.hsp_query_from = 0
        self.hsp_query_to = 0
        self.hsp_hit_from = 0
        self.hsp_hit_to = 0
        self.hsp_identity = 0
        self.hsp_positive = 0
        self.hsp_gaps = 0
        self.hsp_align_len = 0
        self.hsp_qseq = ""
        self.hsp_hseq = ""
        self.hsp_midline = ""


# ==================== DEMONSTRATION CODE ====================

if __name__ == '__main__':
    """
    Demonstration of parser usage.
    Shows how to:
    1. Parse BLAST XML file
    2. Identify queries with no hits (potential pathogenicity factors)
    3. Display results in useful format
    """

    import sys

    # Use command line argument if provided, otherwise use default
    filename = sys.argv[1] if len(sys.argv) > 1 else '/Users/franciscomeersohn/Downloads/blast/bin/MEERSOHNF_blast_out.xml'

    parser = BlastXMLParser()

    # Statistics
    total_queries = 0
    queries_with_hits = 0
    queries_no_hits = 0
    no_hit_queries = []

    print("=" * 80)
    print("BLAST XML Parser Results")
    print("=" * 80)
    print()

    # Parse the file
    for record in parser.parse(filename):
        total_queries += 1

        if record.has_hits:
            queries_with_hits += 1
        else:
            queries_no_hits += 1
            no_hit_queries.append(record)

    # Display summary statistics
    print(f"Total queries processed: {total_queries}")
    print(f"Queries with hits: {queries_with_hits}")
    print(f"Queries with NO hits: {queries_no_hits}")
    print()

    # Display queries with no hits (potential pathogenicity factors!)
    if no_hit_queries:
        print("=" * 80)
        print("QUERIES WITH NO HITS (Potential O157:H7-specific proteins)")
        print("=" * 80)
        print()

        for record in no_hit_queries:
            print(f"Query ID: {record.query_id}")
            print(f"  Definition: {record.query_def}")
            print(f"  Length: {record.query_len} aa")
            print(f"  Message: {record.message if record.message else 'No hits found'}")
            print()

    # Optional: Show example of a query with hits
    print("=" * 80)
    print("SAMPLE OUTPUT: First query with hits")
    print("=" * 80)

    parser2 = BlastXMLParser()
    for record in parser2.parse(filename):
        if record.has_hits:
            print(f"\nQuery: {record.query_id} - {record.query_def}")
            print(f"Length: {record.query_len} aa")
            print(f"Number of hits: {record.hit_count}")
            print("\nTop hit:")
            top_hit = record.hits[0]
            print(f"  Hit ID: {top_hit.hit_id}")
            print(f"  Definition: {top_hit.hit_def}")
            print(f"  Length: {top_hit.hit_len} aa")
            print(f"  Best E-value: {top_hit.best_evalue}")
            if top_hit.hsps:
                best_hsp = top_hit.hsps[0]
                print(f"  Identity: {best_hsp.percent_identity:.1f}%")
                print(f"  Alignment length: {best_hsp.align_len} aa")
            break

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)