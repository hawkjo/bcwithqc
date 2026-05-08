import logging
from Bio import Align
from Bio.Seq import Seq

log = logging.getLogger(__name__)

class CustomBCAligner:
    aligner = Align.PairwiseAligner()
    aligner.wildcard = 'N'
    aligner.mismatch = -1
    aligner.gap_score = -1.1
    aligner.target_left_gap_score = -1.9
    aligner.query_left_gap_score = -1.9
    aligner.target_right_gap_score = 0
        
    def __init__(self, *args, unknown_read_orientation=False):
        """
        Input a list of prefixes that are strings of interest. Ns wild.
        """
        self._unknown_read_orientation = unknown_read_orientation
        self.prefixes = args
        self.full_prefix = ''.join(self.prefixes)
        self.prefix_ends = [sum(len(p) for p in self.prefixes[:i+1]) for i in range(len(self.prefixes))]
        self.prefix_all_Ns = [set(prefix) == set('N') for prefix in self.prefixes]
        self.max_query_len = int(1.5*len(self.full_prefix))
        

    def _log_alignment(self, alignment, seq, orientation):
        """Log the exact Bio.Align pretty-printed alignment for debugging.

        This is guarded by DEBUG level so normal runs do not pay the string
        formatting/log-volume cost unless debugging is explicitly enabled.
        """
        if not log.isEnabledFor(logging.DEBUG):
            return

        log.debug(
            "BC alignment (%s): score=%s, normalized_score=%s, "
            "target_len=%s, query_len=%s, max_query_len=%s\n%s",
            orientation,
            alignment.score,
            alignment.score / len(self.full_prefix),
            len(self.full_prefix),
            len(seq[:self.max_query_len]),
            self.max_query_len,
            format(alignment),
        )

    def find_norm_score_and_key_boundaries(self, seq: Seq):
        """
        Find best alignment and return norm_score=score/alignment_length and boundary positions.
        """
        alignment = self.aligner.align(self.full_prefix, str(seq[:self.max_query_len]))[0]
        orientation = "forward"
        if self._unknown_read_orientation:
            revcompseq = seq.reverse_complement()
            alignment2 = self.aligner.align(self.full_prefix, str(revcompseq[:self.max_query_len]))[0]
            if log.isEnabledFor(logging.DEBUG):
                self._log_alignment(alignment, seq, "forward_candidate")
                self._log_alignment(alignment2, revcompseq, "reverse_complement_candidate")
            if alignment2.score > alignment.score:
                alignment = alignment2
                seq = revcompseq
                orientation = "reverse_complement"

        self._log_alignment(alignment, seq, orientation)

        obs_ends = [None for _ in range(len(self.prefixes))]
        obs_idx = 0
        for i in range(len(alignment.aligned[0])):
            tstart, tend = alignment.aligned[0][i]
            qstart, qend = alignment.aligned[1][i]

            for obs_idx, prefix_end in enumerate(self.prefix_ends[:-1]):
                if obs_ends[obs_idx] is None:
                    # Biopython’s aligned intervals are zero-based and half-open intervals: tstart <= position < tend
                    if tstart <= prefix_end < tend:
                        # log.debug("We are in the 'if tstart <= prefix_end < tend' section")
                        obs_ends[obs_idx] = qstart + prefix_end - tstart

                    elif tstart >= prefix_end:
                        if i == 0 and obs_idx == 0:  # bizarre alignment. discard
                            return None
                        # log.debug("We are in the 'elif tstart >= prefix_end' section")
                        # We have passed or reached the target boundary.
                        # If there was a query insertion between the previous aligned block
                        # and this aligned block, assign that insertion to the block on the left.
                        obs_ends[obs_idx] = qstart

            if obs_ends[-2] is not None:
                break
                    
                
        tstart, tend = alignment.aligned[0][-1]
        qstart, qend = alignment.aligned[1][-1]
        obs_ends[-1] = qend + len(self.full_prefix) - tend
        log.debug("Observed prefix ends at positions: %s", obs_ends)
        return alignment.score/len(self.full_prefix), obs_ends, seq

        
    def find_norm_score_and_pieces(self, seq: Seq, return_seq=False):
        norm_score, obs_ends, seq = self.find_norm_score_and_key_boundaries(seq)
        pieces = [str(seq[:obs_ends[0]])] + [str(seq[obs_ends[i]:obs_ends[i+1]]) for i in range(len(obs_ends)-1)]
        return (norm_score, pieces) if not return_seq else (norm_score, pieces, seq)

    def find_norm_score_pieces_and_boundaries(self, seq: Seq, return_seq=False):
        norm_score, obs_ends, seq = self.find_norm_score_and_key_boundaries(seq)
        pieces = [str(seq[:obs_ends[0]])] + [str(seq[obs_ends[i]:obs_ends[i+1]]) for i in range(len(obs_ends)-1)]
        return (norm_score, pieces, obs_ends) if not return_seq else (norm_score, pieces, obs_ends, seq)
