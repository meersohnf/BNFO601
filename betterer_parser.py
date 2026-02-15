__author__ = 'Paul Fawcett'

from dataclasses import dataclass
from decimal import Decimal
from typing import Tuple
import re

# The two classes below illustrate the modern python approach to creating lightweight data classes
# dataclass is a decorator that auto-creates a class with __init__ and __repr__ methods.
# Type hinting here is especially helpful, although I personally use it sparingly because I'm old and stubborn.
# The structure here implicitly recognizes that for each query sequence, we can have multiple "hits"
# frozen=True here is used to make each record immutable once it is created, which is far safer.

@dataclass(frozen=True)
class Hit:
    subject_id: str
    annotation: str
    evalue: str   # or float, depending on input reliability
    subject_length: int

    # What follows demonstrates two clever tricks every bioinformatician should know
    # When I define a method using the property decorator, I can treat method as if it were a data attribute
    # You don't need to do something like evalue = Hit.evalue_decimal() like when you invoke a method
    # Instead, you can just say evalue_decimal = hit.evalue_decimal just like it was a variable.
    # I am doing this by the way, because some of the evalues reported by blast are too small to record in a Python
    # float.  Doing it this way, I can record both the string representation but also optionally get back a numeric
    # representation using the Decimal type, which has arbitrary precision (but is way slower than floats for math)
    @property
    def evalue_decimal(self): return Decimal(self.evalue)

@dataclass(frozen=True)
class QueryResult:
    query_id: str
    query_annot: str
    query_len: int
    hits: Tuple[Hit, ...]

class BlastParser:

    """
    Parses BLAST (Basic Local Alignment Search Tool) output.

    The BlastParser class is designed to parse BLAST output files and extract relevant
    information such as the query, subject, alignment length, E-value, etc. The parsing
    is facilitated using predefined regular expression patterns. The state-based
    architecture ensures efficient processing of the output in a step-by-step manner
    as the parser transitions between different parsing states. Our parser is literally a state machine.

    Note that nowadays, it would actually be better to use an XML-based parser for this and related tasks.
    BLAST output can be produced directly in XML format using the "-outfmt 5" option.
    As we discovered, BLAST text output is not stable across different versions, but XML output is.
    """

    # Items that are constant across all BlastParser instances, like the regular expressions
    # defining the patterns we are searching for during parsing are well suited to being class attributes.

    # Note that these patterns are "brittle" in the sense that they are very specific to the format of the BLAST output.
    # If the format changes, these patterns will need to be updated accordingly.

    Q_PATTERN = re.compile(r'^Query=\s+(Z\d{4})\s+(.*)')
    S_PATTERN = re.compile(r'^>(\S+)\s+(.*)')
    LEN_PATTERN = re.compile(r'Length=(\d+)')
    EVAL_PATTERN = re.compile(r'Expect\s*=\s*([^,\s]+)')
    NO_HIT_PATTERN = re.compile(r'No hits found')

    def __init__(self):
        """
        Note that in this version, the initializer doesn't do any "work", but rather just sets stuff up.
        This is much cleaner because now a BlastParser instance can be instantiated without
        necessarily having to immediately parse a particular file.
        Style guidelines suggest that initializers should be kept as simple as possible.
        """

        self._reset_query_state()
        self.state = self._scan_for_query  # The initial state of the parser is to scan for the first query
        # Note that we are not invoking the _scan_for_query method here (there are no parenthesis following it)
        # but rather just setting it as an instance variable. self.state is now a reference to a method.

    def parse(self, filename: str):

        """
        This is an important change from BetterParser.py
        Here now we can parse a file using a generator that uses yield to cough up the results
        one record at a time. This allows us to treat the BLAST output file as an iterable...
        not just line-by-line as with a filehandle, but now record-by-record!
        Note that this method is the user interface to the parser. The others following are all private and
        intended for internal use.
        """

        with open(filename) as ifh:  # Here we follow best practices and use a context manager
            for line in ifh:
                result = self.state(line.strip())
                if result is not None:
                    yield result  # "yielding" the result instead of returning sets up the generator.

        # emit the final record if the file ends mid-query
        final = self._finalize_query()
        if final is not None:
            yield final

    # The methods in the following block all correspond to an internal state of a parser when using
    # our state-based parsing logic. These will be selected as appropriate based on the current state of the parser
    # as it transitions between states.
    # Note here a new architecture! Now a state is a function that always consumes a line and optionally emits a record.
    # We are better off explicitly passing lines rather than relying on an instance variable.
    # Note also that we now use just a single underscore for method privacy, reflecting modern python best practices
    # Now the double underscore is used mostly just to protect the method from accidental overwriting of the method
    # when something similarly named is used in a subclass (a so-called "subclass collision").
    # So it is now only commonly used when a class is intended to be subclassed using python's inheritance mechanism.

    def _scan_for_query(self, line: str):
        m = self.Q_PATTERN.search(line)
        if m:
            result = self._finalize_query()
            self.query_id, self.query_annot = m.groups()
            self.state = self._extend_query # Change state to extend this query
            return result
        return None # Making this explicit for clarity

    def _extend_query(self, line: str):
        m = self.LEN_PATTERN.search(line)
        if m:
            self.query_len = int(m.group(1))
            self.state = self._scan_for_subject  # Change state to scan for subjects
        else:
            self.query_annot += ' ' + line.strip()

    def _scan_for_subject(self, line: str):
        if self.NO_HIT_PATTERN.search(line):
            self.hits = []
            result = self._finalize_query()
            self.state = self._scan_for_query
            return result

        m = self.S_PATTERN.search(line)
        if m:
            self.subject_id, self.subject_annot = m.groups()
            self.state = self._extend_subject
            return None

        return self._scan_for_query(line)

    def _extend_subject(self, line: str):
        m = self.LEN_PATTERN.search(line)
        if m:
            self.subject_length = int(m.group(1))
            self.state = self._scan_for_eval # We found the Length pattern, subject is complete
        elif line.strip():  # This is just going to skipp blank lines altogether
            self.subject_annot += ' ' + line.strip()  # otherwise we must be extending the subject

    def _scan_for_eval(self, line: str):
        m = self.EVAL_PATTERN.search(line)
        if m:
            hit = Hit(
                subject_id=self.subject_id,
                annotation=self.subject_annot,
                evalue=m.group(1),
                subject_length=self.subject_length
            )
            self.hits.append(hit)
            self.state = self._scan_for_subject

    #  The following two methods are really just helpers and don't correspond to an internal state of the parser.

    def _finalize_query(self):
        if not self.query_id:
            return None

        result = QueryResult(
            query_id=self.query_id,
            query_annot=self.query_annot,
            query_len=self.query_len,
            hits=tuple(self.hits)  # The list is converted to a tuple to make it immutable once everything is found
            # This guarantees against the silent corruption of previously emitted records and helps guarantee safe code
        )
        self._reset_query_state()
        return result

    def _reset_query_state(self):
        self.query_id = ''
        self.query_annot = ''
        self.query_len = 0
        self.subject_length = 0
        self.subject_id = ''
        self.subject_annot = ''
        self.hits = []

# Here's an example of how to use the new parser.
# Another improvement is that the new parser only produces records and never prints the results!
# That division of responsibilities is much better. The parser does the actual parsing, and the caller is responsible
# for whatever they want to do with the results.

if __name__ == '__main__':  # This usage example will only run if this file is executed directly, but not if imported

    parser = BlastParser()
    print('Here are the queries that had no hits:')
    for record in parser.parse('K12_edl_results.txt'):  # Note how calling parser.parse can be treated as an iterable
        if record.hits:  # Only interested in queries with no hits
            print(record)
