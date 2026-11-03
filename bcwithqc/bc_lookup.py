from collections import defaultdict


class bc_lookup:
    """
    Chunk-based barcode candidate lookup.

    The goal is to quickly reduce a query barcode to a small set of possible
    reference barcodes by matching exact chunks instead of scanning the full
    whitelist every time.
    """

    def __init__(self, reference_barcodes, allowed_errors):

        if not reference_barcodes:
            raise ValueError("reference_barcodes must not be empty.")

        if not isinstance(allowed_errors, int) or allowed_errors < 0:
            raise ValueError("allowed_errors must be a non-negative integer.")

        self.reference_barcodes = set(reference_barcodes)
        self.allowed_errors = allowed_errors

        # Get unique lengths of reference barcodes
        self.possible_barcode_lengths = {len(bc) for bc in reference_barcodes}

        self.dict_of_ref_break_functions = self._build_dict_of_ref_break_functions()
        self.query_break_function = self._build_query_break_function()

        # list of dictionaries:
        # [
        #   {chunk_seq: {barcode1, barcode2, ...}},
        #   {chunk_seq: {barcode1, barcode2, ...}},
        #   ...
        # ]
        self.list_of_lookup_dicts = self._build_list_of_lookup_dicts(reference_barcodes, )




    def _build_dict_of_ref_break_functions(self):
        """
        Build a dictionary of functions that break a barcode of the given length into chunks.

        keys = barcode length
        values = function that takes a barcode and returns chunks
        """

        dict_of_ref_break_functions = {}

        for barcode_length in self.possible_barcode_lengths:

            kmer = barcode_length // (self.allowed_errors + 1) 
            leftover_bases = barcode_length % (self.allowed_errors + 1)
            n_pieces = self.allowed_errors + 1

            def ref_break(barcode, kmer=kmer, leftover_bases=leftover_bases, n_pieces=n_pieces):

                chunks = []
                start = 0

                for piece_i in range(n_pieces):
                    chunk_length = kmer

                    if piece_i < leftover_bases:
                        chunk_length += 1

                    end = start + chunk_length
                    chunks.append(barcode[start:end])
                    start = end

                return tuple(chunks)

            dict_of_ref_break_functions[barcode_length] = ref_break

        return dict_of_ref_break_functions


    def _make_query_boundaries(self, query_len, ref_len, allowed_errors):
        """
        Precompute all query slices that could correspond to a query length
        and reference barcode length.

        Returns
        -------
        tuple of:
            (piece_i, ((start1, end1), (start2, end2), ...))
        """
        n_pieces = allowed_errors + 1

        kmer, leftover_bases = divmod(ref_len, n_pieces)

        boundaries_by_piece = {}

        ref_start = 0

        # Get reference pieces and add leftover bases to the first pieces as needed.
        for piece_i in range(n_pieces):
            if piece_i < leftover_bases:
                piece_len = kmer + 1
            else:
                piece_len = kmer

            ref_end = ref_start + piece_len

            breakpoints = set()

            # Shift the query slice left and right by allowed_errors to account for possible insertions/deletions.
            for shift in range(-allowed_errors, allowed_errors + 1):
                query_start = ref_start + shift
                query_end = query_start + piece_len

                # Discard slices that fall outside the observed query barcode.
                if query_start < 0:
                    continue
                # Discard slices that no longer match reference piece length.
                if query_end > query_len:
                    continue

                breakpoints.add((query_start, query_end))

            if breakpoints:
                boundaries_by_piece[piece_i] = tuple(sorted(breakpoints))

            ref_start = ref_end

        # Turn mutable dict into a sorted tuple of tuples for better itteration. 
        # before:
        # boundaries_by_piece = {
        #     0: {(0, 7), (1, 8)},
        #     1: {(7, 14), (8, 15)},
        #     2: {(14, 20)}
        # }
        # After:
        #  (
        #     (0, ((0, 7), (1, 8))),
        #     (1, ((7, 14), (8, 15))),
        #     (2, ((14, 20),)),
        # )
        return tuple(
            (piece_i, boundaries_by_piece[piece_i])
            for piece_i in sorted(boundaries_by_piece)
        )

    def _build_query_break_function(self):
        allowed_errors = self.allowed_errors
        possible_lengths = tuple(sorted(self.possible_barcode_lengths))

        min_query_len = max(min(possible_lengths) - allowed_errors, 1)
        max_query_len = max(possible_lengths) + allowed_errors

        breaks_by_query_len = {}

        for query_len in range(min_query_len, max_query_len + 1):
            # Temporary build-time structure:
            # piece_i -> set of unique (start, end) breakpoints
            breakpoints_by_piece = {}

            for ref_len in possible_lengths:
                if abs(query_len - ref_len) > allowed_errors:
                    continue

                boundaries_by_piece = self._make_query_boundaries(
                    query_len=query_len,
                    ref_len=ref_len,
                    allowed_errors=allowed_errors,
                )

                for piece_i, breakpoints in boundaries_by_piece:
                    # Initialize the set for this piece if it doesn't exist yet
                    if piece_i not in breakpoints_by_piece:
                        breakpoints_by_piece[piece_i] = set()
                    # Update the set of breakpoints for this piece with the new breakpoints
                    breakpoints_by_piece[piece_i].update(breakpoints)

            if breakpoints_by_piece:
                # Freeze into a simple runtime structure:
                # query_len -> ((piece_i, ((start, end), ...)), ...)
                breaks_by_query_len[query_len] = tuple(
                    (piece_i, tuple(sorted(breakpoints)))
                    for piece_i, breakpoints in sorted(breakpoints_by_piece.items())
                )

        def query_break(barcode, breaks_by_query_len=breaks_by_query_len):
            piece_plans = breaks_by_query_len.get(len(barcode))
            if piece_plans is None:
                return ()

            out = []

            for piece_i, breakpoints in piece_plans:
                kmers = set()

                for start, end in breakpoints:
                    kmers.add(barcode[start:end])

                out.append((piece_i, kmers))

            return tuple(out)

        return query_break

    def _build_list_of_lookup_dicts(self, reference_barcodes):
        """
        Build list of lookup dictionaries for each kmer number.

        Each dictionary has keys = chunk sequence, values = set of reference barcodes
        that contain that chunk in the correct position.
        """

        list_of_lookup_dicts = [defaultdict(set) for _ in range(self.allowed_errors + 1)]

        for barcode in reference_barcodes:
            barcode_length = len(barcode)

            ref_break_func = self.dict_of_ref_break_functions[barcode_length]
            chunks = ref_break_func(barcode)

            for kmer_number, chunk in enumerate(chunks):
                list_of_lookup_dicts[kmer_number][chunk].add(barcode)

        return list_of_lookup_dicts

    def get_candidate_barcodes(self, query_barcode):
        """
        Return possible reference barcodes for a query barcode.
        """
        candidate_barcodes = set()

        query_kmer_chunk_sets = self.query_break_function(query_barcode)

        lookup_dicts = self.list_of_lookup_dicts
        n_lookup_dicts = len(lookup_dicts)

        for piece_i, query_chunks in query_kmer_chunk_sets:
            if piece_i >= n_lookup_dicts:
                # This should never happen if query_break_function was built correctly.
                raise ValueError(
                    f"Query barcode {query_barcode} produced piece_i={piece_i}, "
                    f"but only {n_lookup_dicts} lookup dictionaries exist."
                )

            lookup_dict = lookup_dicts[piece_i]

            for query_chunk in query_chunks:
                candidate_barcodes.update(
                    lookup_dict.get(query_chunk, ())
                )

        return candidate_barcodes