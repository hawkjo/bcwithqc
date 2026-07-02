# tests/test_14_bc_lookup.py

import os
import random

import pytest

from bcwithqc.bc_lookup import bc_lookup


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def small_lookup():
    return bc_lookup(
        reference_barcodes=[
            "AAACCC",
            "AAAGGG",
            "TTTCCC",
            "GGGAAA",
        ],
        allowed_errors=2,
    )


def test_ref_break_function_dict_has_expected_lengths(small_lookup):
    assert small_lookup.possible_barcode_lengths == {6}
    assert set(small_lookup.dict_of_ref_break_functions.keys()) == {6}


def test_ref_break_function_splits_length_6_barcode_correctly(small_lookup):
    break_func = small_lookup.dict_of_ref_break_functions[6]

    assert break_func("AAACCC") == ("AA", "AC", "CC")
    assert break_func("AAAGGG") == ("AA", "AG", "GG")
    assert break_func("TTTCCC") == ("TT", "TC", "CC")
    assert break_func("GGGAAA") == ("GG", "GA", "AA")


def test_lookup_dicts_are_built_correctly(small_lookup):
    """
    Check exact contents of the lookup dictionaries.

    With allowed_errors = 2, barcodes are split into 3 chunks.
    Therefore list_of_lookup_dicts has 3 dictionaries:
        index 0 = first chunk
        index 1 = second chunk
        index 2 = third chunk
    """

    assert len(small_lookup.list_of_lookup_dicts) == 3

    assert small_lookup.list_of_lookup_dicts[0]["AA"] == {
        "AAACCC",
        "AAAGGG",
    }
    assert small_lookup.list_of_lookup_dicts[0]["TT"] == {
        "TTTCCC",
    }
    assert small_lookup.list_of_lookup_dicts[0]["GG"] == {
        "GGGAAA",
    }

    assert small_lookup.list_of_lookup_dicts[1]["AC"] == {
        "AAACCC",
    }
    assert small_lookup.list_of_lookup_dicts[1]["AG"] == {
        "AAAGGG",
    }
    assert small_lookup.list_of_lookup_dicts[1]["TC"] == {
        "TTTCCC",
    }
    assert small_lookup.list_of_lookup_dicts[1]["GA"] == {
        "GGGAAA",
    }

    assert small_lookup.list_of_lookup_dicts[2]["CC"] == {
        "AAACCC",
        "TTTCCC",
    }
    assert small_lookup.list_of_lookup_dicts[2]["GG"] == {
        "AAAGGG",
    }
    assert small_lookup.list_of_lookup_dicts[2]["AA"] == {
        "GGGAAA",
    }


def test_query_break_function_returns_expected_shifted_chunks(small_lookup):
    query_chunks = small_lookup.query_break_function("AAACCC")

    assert dict(query_chunks) == {
        0: {"AA", "AC"},
        1: {"AA", "AC", "CC"},
        2: {"AC", "CC"},
    }

def test_get_candidate_barcodes_returns_exact_expected_candidates(small_lookup):
    """
    Query AAACCC should return candidates that share at least one allowed
    shifted chunk in the correct kmer position.
    """

    candidates = small_lookup.get_candidate_barcodes("AAACCC")

    assert candidates == {
        "AAACCC",
        "AAAGGG",
        "TTTCCC",
    }


def test_get_candidate_barcodes_can_return_single_exact_candidate(small_lookup):
    candidates = small_lookup.get_candidate_barcodes("GGGAAA")

    assert candidates == {
        "GGGAAA",
    }


def test_random_query_candidate_percentage_is_written_to_file():
    """
    Load a larger reference barcode list and estimate how many reference
    barcodes are returned on average for random query sequences.

    This writes a TSV table with one row per allowed_errors value.
    """

    ref_barcodes_fpath = os.path.join(
        SCRIPT_DIR,
        "../tests/class_bc_lookup/barcodes.txt",
    )

    assert os.path.exists(ref_barcodes_fpath), (
        f"Reference barcode file does not exist: {ref_barcodes_fpath}"
    )

    with open(ref_barcodes_fpath) as handle:
        ref_barcodes = [
            line.strip()
            for line in handle
            if line.strip()
        ]

    assert len(ref_barcodes) > 0

    rng = random.Random(1)
    nucleotides = ("A", "C", "G", "T")

    n_queries = 1000
    random_query_barcodes = []
    query_lengths = []

    for _ in range(n_queries):
        # ca. 70% length 20, remaining 30% spread over 18, 19, 21, 22
        if rng.random() < 0.70:
            query_length = 20
        else:
            query_length = rng.choice([18, 19, 21, 22])

        query_lengths.append(query_length)

        query_barcode = "".join(
            rng.choice(nucleotides)
            for _ in range(query_length)
        )

        random_query_barcodes.append(query_barcode)

    total_ref_barcodes = len(set(ref_barcodes))

    results = []

    for allowed_errors in range(0, 11):
        lookup = bc_lookup(
            reference_barcodes=ref_barcodes,
            allowed_errors=allowed_errors,
        )

        candidate_percentages = []

        for query_barcode in random_query_barcodes:
            candidates = lookup.get_candidate_barcodes(query_barcode)

            candidate_percentage = (
                len(candidates) / total_ref_barcodes
            ) * 100

            candidate_percentages.append(candidate_percentage)

        mean_candidate_percentage = (
            sum(candidate_percentages) / len(candidate_percentages)
        )

        results.append({
            "allowed_errors": allowed_errors,
            "n_reference_barcodes": total_ref_barcodes,
            "n_random_queries": n_queries,
            "mean_candidate_percentage": mean_candidate_percentage,
            "min_candidate_percentage": min(candidate_percentages),
            "max_candidate_percentage": max(candidate_percentages),
        })

    output_fpath = os.path.join(
        os.path.dirname(ref_barcodes_fpath),
        "bc_lookup_random_query_candidate_percentage_by_allowed_errors.tsv",
    )

    query_length_counts = {
        query_length: query_lengths.count(query_length)
        for query_length in sorted(set(query_lengths))
    }

    with open(output_fpath, "w") as handle:
        handle.write("# bc_lookup random query candidate percentage summary\n")
        handle.write(f"# reference_barcode_file\t{ref_barcodes_fpath}\n")
        handle.write(f"# n_reference_barcodes\t{total_ref_barcodes}\n")
        handle.write(f"# n_random_queries\t{n_queries}\n")

        for query_length, count in query_length_counts.items():
            handle.write(f"# query_length_{query_length}_count\t{count}\n")

        handle.write("\n")

        handle.write(
            "allowed_errors\t"
            "k\t"
            "n_reference_barcodes\t"
            "n_random_queries\t"
            "mean_candidate_percentage\t"
            "min_candidate_percentage\t"
            "max_candidate_percentage\n"
        )

        for result in results:
            handle.write(
                f"{result['allowed_errors']}\t"
                f"{20//(result['allowed_errors'] + 1)}\t"
                f"{result['n_reference_barcodes']}\t"
                f"{result['n_random_queries']}\t"
                f"{result['mean_candidate_percentage']:.6f}\t"
                f"{result['min_candidate_percentage']:.6f}\t"
                f"{result['max_candidate_percentage']:.6f}\n"
            )

    assert os.path.exists(output_fpath)

    for result in results:
        assert 0 <= result["mean_candidate_percentage"] <= 100
        assert 0 <= result["min_candidate_percentage"] <= 100
        assert 0 <= result["max_candidate_percentage"] <= 100