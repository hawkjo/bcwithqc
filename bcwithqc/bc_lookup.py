from collections import defaultdict


class bc_lookup:
    """
    Chunk-based barcode candidate lookup.

    The goal is to quickly reduce a query barcode to a small set of possible
    reference barcodes by matching exact chunks instead of scanning the full
    whitelist every time.
    """

    def __init__(self, reference_barcodes, allowed_errors):
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

            def ref_break(barcode, barcode_length=barcode_length):
                kmer = barcode_length // (self.allowed_errors + 1)
                leftover_bases = barcode_length % (self.allowed_errors + 1)

                chunks = []
                start = 0

                for _ in range(self.allowed_errors + 1):
                    chunk_length = kmer

                    if leftover_bases > 0:
                        chunk_length += 1
                        leftover_bases -= 1

                    end = start + chunk_length
                    chunks.append(barcode[start:end])
                    start = end

                return tuple(chunks)

            dict_of_ref_break_functions[barcode_length] = ref_break

        return dict_of_ref_break_functions

    def _build_query_break_function(self):
        """
        Build a function that breaks a query barcode into possible shifted chunks.

        The returned function takes one query barcode and returns a tuple of sets:
            index 0 = possible chunks for kmer position 0
            index 1 = possible chunks for kmer position 1
            etc.
        """

        def query_break(barcode):

            # Create list of lists:
            # index 0 = chunks for kmer 0
            # index 1 = chunks for kmer 1
            # etc.
            kmer_chunks = [set() for _ in range(self.allowed_errors + 1)]

            for barcode_length in self.possible_barcode_lengths:
                kmer = barcode_length // (self.allowed_errors + 1)
                leftover_bases = barcode_length % (self.allowed_errors + 1)

                start = 0

                for kmer_number in range(self.allowed_errors + 1):
                    chunk_length = kmer

                    if leftover_bases > 0:
                        chunk_length += 1
                        leftover_bases -= 1

                    end = start + chunk_length

                    # Shift window left and right
                    for shift in range(-self.allowed_errors, self.allowed_errors + 1):
                        shifted_start = start + shift
                        shifted_end = end + shift

                        # Only keep full-length chunks inside barcode boundaries
                        if shifted_start < 0:
                            continue

                        if shifted_end > len(barcode):
                            continue

                        shifted_chunk = barcode[shifted_start:shifted_end]

                        # Extra safety check
                        if len(shifted_chunk) == chunk_length:
                            kmer_chunks[kmer_number].add(shifted_chunk)

                    start = end

            return tuple(kmer_chunks)

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

        for kmer_number, query_chunks in enumerate(query_kmer_chunk_sets):
            if kmer_number < len(self.list_of_lookup_dicts):

                for query_chunk in query_chunks:
                    candidate_barcodes.update(
                        self.list_of_lookup_dicts[kmer_number].get(query_chunk, set())
                    )

            else:
                # This should never happen.
                raise ValueError(
                    f"Query barcode {query_barcode} produced a higher kmer_number "
                    f"than should be possible."
                )

        return candidate_barcodes