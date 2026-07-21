import freebarcodes.decode
import logging
import sys

from .misc import DistanceThresh
from .bc_lookup import bc_lookup
from collections import Counter

log = logging.getLogger(__name__)

class BCDecoder:
    def __init__(self, bc_whitelist, bc_maxdist):
        self.bcs = bc_whitelist
        self.bcs_set = set(self.bcs)
        self.bc_maxdist = bc_maxdist
        self.bc_len = len(self.bcs[0])
        self._distfun = DistanceThresh("levenshtein", bc_maxdist)

        # I removed this assertion because it should no longer be required, as we can now handle barcodes of different lengths.
        # assert all(len(bc) == self.bc_len for bc in self.bcs)

        self.k = min(map(len, self.bcs)) // (self.bc_maxdist + 1)
        if self.k > 2:
            self.bc_lookup = bc_lookup(
                reference_barcodes=self.bcs,
                allowed_errors=self.bc_maxdist,
            )

    def _candidate_bcs(self, raw_bc):
        """
        Return candidate whitelist barcodes for a raw barcode.

        The lookup only prefilters candidates. Final validity is still decided
        by DistanceThresh.
        """
        return self.bc_lookup.get_candidate_barcodes(raw_bc)

    def decode(self, raw_bc):
        if raw_bc in self.bcs_set:
            return raw_bc
        if self.k > 2:
            candidates = self._candidate_bcs(raw_bc)
            if not candidates:
                return None
            dists_and_scores = [(dist, bc) for bc in candidates if (dist := self._distfun(raw_bc, bc)) is not False]

            if not dists_and_scores:
                return None

            min_dist, bc = min(dists_and_scores)
            if sum(dist == min_dist for dist, _ in dists_and_scores) == 1:
                return bc

            return None
        else:
            dists_and_scores = [(dist, bc) for bc in self.bcs if (dist := self._distfun(raw_bc, bc)) is not False]
            if not len(dists_and_scores):
                return None

            min_dist, bc = min(dists_and_scores)

            if sum(dist == min_dist for dist, _ in dists_and_scores) > 1:
                return None

            return bc

    def decode_with_status(self, raw_bc):
        if raw_bc in self.bcs_set:
            return raw_bc, "exact", None

        if self.k > 2:
            candidates = self._candidate_bcs(raw_bc)
            if not candidates:
                return None, "no_match", None
            
            dists_and_scores = [(dist, bc) for bc in candidates if (dist := self._distfun(raw_bc, bc)) is not False]

            if not dists_and_scores:
                return None, "no_match", None

            min_dist = min(dist for dist, _ in dists_and_scores)
            min_dist_bcs = [bc for dist, bc in dists_and_scores if dist == min_dist]

            if len(min_dist_bcs) > 1:
                return None, "ambiguous", min_dist_bcs

            return min_dist_bcs[0], "corrected", None
        else:
            dists_and_scores = [(dist, bc) for bc in self.bcs if (dist := self._distfun(raw_bc, bc)) is not False]

            if not dists_and_scores:
                return None, "no_match", None

            min_dist = min(dist for dist, _ in dists_and_scores)
            min_dist_bcs = [bc for dist, bc in dists_and_scores if dist == min_dist]

            if len(min_dist_bcs) > 1:
                return None, "ambiguous", min_dist_bcs

            return min_dist_bcs[0], "corrected", None

class SBCDecoder:
    def __init__(self, sbc_whitelist, sbc_maxdist, sbc_reject_delta):
        self.sbcs = sbc_whitelist
        # This will fail for whitelist barcodes of different lengths, might need to fix that. 
        # This branch does currently not support barcodes of different lengths!
        self.sbc_len = self._validate_sbc_lengths()
        
        self.sbc_maxdist = sbc_maxdist
        self.sbc_reject_delta = sbc_reject_delta
        self.sbcd = freebarcodes.decode.FreeDivBarcodeDecoder()
        self.sbcd.build_codebook_from_random_codewords(self.sbcs, self.sbc_maxdist, self.sbc_reject_delta)

    def decode(self, raw_sbc):
        sbc = self.sbcd.decode(raw_sbc)
        return sbc if isinstance(sbc, str) else None

    def decode_with_status(self, raw_sbc):
        result = self.sbcd.decode(raw_sbc)

        if isinstance(result, str):
            if result == raw_sbc:
                return result, "exact", None
            return result, "corrected", None

        if result is None:
            return None, "no_match", None

        if isinstance(result, int) and result < 0:
            conflict_level = -result
            # Call find_conflicting_sbcs, only if there is a conflict. 
            conflicting_sbcs = self.find_conflicting_sbcs(raw_sbc, conflict_level)
            return None, "ambiguous", conflicting_sbcs

        raise ValueError(f"Unexpected decoder result for {raw_sbc!r}: {result!r}")
    
    def find_conflicting_sbcs(self, raw_sbc, conflict_level):
        """
        Return whitelist SBCs that plausibly participate in the ambiguous conflict.

        conflict_level is the positive value corresponding to the negative integer
        returned by FreeDivBarcodeDecoder.decode().
        """
        # calculate conflict radius, mirroring the decoder build 
        # We only care about whitelist barcodes within this free divergence
        max_conflict_radius = conflict_level + self.sbc_reject_delta

        # Calculate the free divergence to every whitelist barcode
        # WARNING: This might reduce the speed at which it runs slightly,
        # but since it only triggers on conflicts, it should not be too bad. 
        dists_and_sbcs = [
            (freebarcodes.editmeasures.free_divergence(raw_sbc, sbc), sbc)
            for sbc in self.sbcs
        ]

        # Identify all barcodes within the max_conflic_radius
        conflicting = [
            (dist, sbc)
            for dist, sbc in dists_and_sbcs
            if dist <= max_conflict_radius
        ]
        # Sort first by distance, then alphabetically by barcode.
        conflicting.sort(key=lambda x: (x[0], x[1]))

        # Here we throw away the "distance" and keep only the conflicting barcodes
        # Might be interesting to keep the distance too. 
        return [sbc for dist, sbc in conflicting]

    def _validate_sbc_lengths(self):
        """
        Validate that all whitelist SBCs have the same length.

        SBCDecoder currently requires fixed-length whitelist barcodes. If
        different lengths are present, log the length distribution and every
        barcode that does not have the primary length before raising an error.

        Returns
        -------
        int
            The common SBC length when validation succeeds.
        """
        if not self.sbcs:
            raise ValueError(
                "SBCDecoder received an empty whitelist. "
                "At least one barcode sequence is required."
            )

        length_counts = Counter(len(sbc) for sbc in self.sbcs)

        # The primary length is the most frequently occurring barcode length.
        # Counter.most_common() resolves ties using first occurrence order.
        primary_length, primary_count = length_counts.most_common(1)[0]

        # If all bcs are the same length, we return the length and move on. 
        if len(length_counts) == 1:
            return primary_length

        # Now we are in the failstate, with bcs of different lengths. 
        total_count = len(self.sbcs)

        log.error(
            "SBCDecoder (Free Divergence) received whitelist barcodes of different lengths. ",
            "Barcodes of different lengths are currently only supported by BCDecoder (levensthein), not SBCDecoder."
        )
        log.error(
            "Barcode length distribution (%d barcodes total):",
            total_count,
        )

        for length, count in sorted(length_counts.items()):
            log.error(
                "  length %d: %d barcode(s)",
                length,
                count,
            )

        log.error(
            "Primary barcode length is %d (%d of %d barcodes).",
            primary_length,
            primary_count,
            total_count,
        )
        log.error(
            "Barcodes that do not have the primary length:"
        )

        incompatible_count = 0

        for index, sbc in enumerate(self.sbcs):
            sbc_length = len(sbc)

            if sbc_length != primary_length:
                incompatible_count += 1
                log.error(
                    "  index=%d, length=%d, sequence=%r",
                    index,
                    sbc_length,
                    sbc,
                )

        raise ValueError(
            "SBCDecoder (Free Divergence) requires all whitelist barcodes to have the same "
            f"length. The primary length is {primary_length}, but "
            f"{incompatible_count} of {total_count} barcodes have a different "
            "length. See the preceding log messages for the length distribution "
            "and incompatible barcode sequences. Variable-length whitelist "
            "barcodes are currently supported by BCDecoder (levensthein), but not SBCDecoder."
        )

