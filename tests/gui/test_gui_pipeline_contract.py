"""Tier F: bcwithqc itself agrees with what the GUI writes and describes.

Opt-in: run with `pytest tests/gui -m bcwithqc`. bcwithqc must be importable
in the same environment (pip install -e ".[test]"). The tests run bcwithqc's
CLI in a subprocess on reads from its own simulate_reads, so no STAR is
needed: preprocess --output-format-bam writes a stand-in BAM that count
accepts. The GUI side runs in Chromium only.

Cases about current bcwithqc limitations are xfail(strict=True). When
bcwithqc is fixed they pass unexpectedly, which fails the run and says what
to change in the GUI (usually: delete a temporary warning).
"""

import copy
import itertools
import json
import os
import random
import subprocess
import sys

import pytest
from gui_support import (
    EXAMPLES_DIR,
    SAMPLE_BARCODES,
    SAMPLE_CONSTANT,
    barcode_list,
    constant_region,
    layout_blocks,
    make_config,
    random_block,
    read_config,
    without_names,
)

pytest.importorskip("bcwithqc")
pytest.importorskip("pysam")

pytestmark = [pytest.mark.bcwithqc, pytest.mark.only_browser("chromium")]

SAMPLE_BARCODES_8 = ["TCGCCTTA", "CTAGTACG", "TTCTGCCT", "GCTCAGGA"]


# ---------------------------------------------------------------------------
# Running bcwithqc


@pytest.fixture(scope="session")
def bcwithqc_env(tmp_path_factory):
    """Environment for bcwithqc subprocesses: no bytecode files, matplotlib's cache in a temporary folder."""
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["MPLCONFIGDIR"] = str(tmp_path_factory.mktemp("matplotlib"))
    return env


def run_bcwithqc(args, folder, env):
    return subprocess.run(
        [sys.executable, "-m", "bcwithqc", *args],
        check=False,
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        timeout=900,
    )


def write_json(config, path):
    path.write_text(json.dumps(config, indent=4) + "\n", encoding="utf-8")
    return path


def zero_errors(config):
    """A copy with maxerrors 0 everywhere, so simulate_reads adds no errors."""
    exact = copy.deepcopy(config)
    for key in ("barcode_struct_r1", "barcode_struct_r2"):
        for block in exact.get(key, {}).get("blocks", []):
            block["maxerrors"] = 0
    return exact


def simulate(config, folder, env, nreads, seed=1, extra=()):
    """Simulate error-free reads for a config; return the folder holding only the FASTQs."""
    config_path = write_json(zero_errors(config), folder / "simulate_config.json")
    reads = folder / "reads"
    reads.mkdir()
    result = run_bcwithqc(
        [
            "simulate_reads",
            f"--config={config_path}",
            f"--output-dir={reads}",
            f"--nreads={nreads}",
            "--error-probability=-1",
            f"--seed={seed}",
            *extra,
        ],
        folder,
        env,
    )
    assert result.returncode == 0, result.stderr
    for log in reads.glob("*.log"):
        log.rename(folder / log.name)
    return reads


def preprocess(config_text, folder, reads, env):
    config_path = folder / "config.json"
    config_path.write_text(config_text, encoding="utf-8")
    output = folder / "pre"
    result = run_bcwithqc(
        [
            "preprocess",
            str(reads),
            f"--config={config_path}",
            f"--output-dir={output}",
            "--threads=1",
            "--output-format-bam",
        ],
        folder,
        env,
    )
    return output, result


def count(folder, output, env):
    return run_bcwithqc(
        [
            "count",
            str(output),
            f"--STAR-output-dir={output}",
            f"--config={folder / 'config.json'}",
            f"--output-dir={output}",
            "--threads=1",
            "--keep-intermediary",
        ],
        folder,
        env,
    )


def read_fastq(path):
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    records = []
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return records
            sequence = handle.readline().rstrip("\n")
            handle.readline()
            handle.readline()
            records.append((header[1:].split()[0], sequence))


def bam_tags(path):
    import pysam

    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return [
            (read.query_name, dict(read.get_tags()))
            for read in bam.fetch(until_eof=True)
        ]


def truth(name, blocks):
    """The true value of each block and the tail length, from a simulate_reads read name."""
    tokens = name.split("_")
    values = {}
    position = 1
    for block in blocks:
        values[block["blockname"]] = tokens[position]
        position += 2 if block["blocktype"] in ("barcodeList", "constantRegion") else 1
    tail = tokens[position]
    assert tail.startswith("randTail"), name
    return values, int(tail[len("randTail") :])


# ---------------------------------------------------------------------------
# The GUI side


def gui_export(core, config):
    """The text the GUI writes for a config, and its description of bcwithqc's outputs."""
    return core.run(
        """(core, arg) => {
            const model = core.config_from_json(arg).model;
            return {text: core.config_to_json(model), outputs: core.describe_outputs(model), issues: core.validate_config(model)};
        }""",
        config,
    )


def expected_tags(outputs, values_by_read):
    """CB, CR, UR and FL as describe_outputs says bcwithqc builds them, from the true block values."""

    def joined(tag):
        per_read = []
        for position, read_index in enumerate(outputs["reads_with_blocks"]):
            names = outputs[tag]["pieces"][position]
            per_read.append(
                ".".join(values_by_read[read_index][name] for name in names)
            )
        return ".".join(per_read)

    fl = sum(
        len(values_by_read[read_index][name])
        for position, read_index in enumerate(outputs["reads_with_blocks"])
        for name in outputs["FL"]["blocks"][position]
    )
    return {"CB": joined("CB"), "CR": joined("CR"), "UR": joined("UR"), "FL": fl}


def expected_alignment_sequence(file_info, sequence, values, tail):
    if file_info["kind"] == "whole":
        return sequence
    pieces = [values[name] for name in file_info["parts"]]
    if file_info["rest"]:
        pieces.append(sequence[len(sequence) - tail :] if tail else "")
    return "".join(pieces)


def fastq_files(reads):
    files = sorted(reads.glob("*.txt.gz"))
    assert files, f"no simulated reads in {reads}"
    return files


def run_end_to_end(core, config, folder, env, nreads=400):
    """Export a config through the GUI, run bcwithqc on error-free reads, and compare every read."""
    exported = gui_export(core, config)
    assert [issue for issue in exported["issues"] if issue["level"] == "error"] == []
    outputs = exported["outputs"]
    gui_config = json.loads(exported["text"])
    read_keys = [
        key for key in ("barcode_struct_r1", "barcode_struct_r2") if key in gui_config
    ]
    blocks_by_read = {
        index: gui_config[key]["blocks"] for index, key in enumerate(read_keys)
    }

    reads = simulate(gui_config, folder, env, nreads)
    simulated = {}
    for read_index, path in enumerate(fastq_files(reads)):
        for name, sequence in read_fastq(path):
            values, tail = truth(name, blocks_by_read[read_index])
            simulated.setdefault(name.split("_")[0], {})[read_index] = (
                name,
                sequence,
                values,
                tail,
            )

    output, result = preprocess(exported["text"], folder, reads, env)
    assert result.returncode == 0, result.stderr[-3000:]
    for file_info in outputs["alignment_files"]:
        if file_info["kind"] == "none":
            continue  # bcwithqc still writes empty records (notes #empty-alignment-files)
        suffix = f"_{file_info['read'] + 1}.fq" if len(read_keys) == 2 else ".fq"
        (path,) = [p for p in output.glob("sans_bc_*.fq") if p.name.endswith(suffix)]
        records = read_fastq(path)
        assert len(records) == nreads
        for name, sequence in records:
            _, raw, values, tail = simulated[name.split("_")[0]][file_info["read"]]
            assert sequence == expected_alignment_sequence(
                file_info, raw, values, tail
            ), name

    result = count(folder, output, env)
    assert result.returncode == 0, result.stderr[-3000:]
    records = bam_tags(output / "with_bc_umi.sorted.bam")
    seen = set()
    for name, tags in records:
        index = name.split("_")[0]
        seen.add(index)
        values_by_read = {
            read_index: item[2] for read_index, item in simulated[index].items()
        }
        want = expected_tags(outputs, values_by_read)
        if not outputs["UR"]["written"]:
            del want[
                "UR"
            ]  # bcwithqc writes it empty today: test_no_umi_block_writes_no_ur_tag
        assert {tag: tags.get(tag) for tag in want} == want, name
        assert ("UB" in tags) is outputs["UB"]["written"], name
    assert seen == set(simulated), "every simulated read reaches the BAM"


# ---------------------------------------------------------------------------
# F1: decoder agreement


def bcwithqc_decoder_labels(blocks):
    from bcwithqc.bc_decoders import SBCDecoder
    from bcwithqc.misc import build_bc_decoders

    decoders = build_bc_decoders(copy.deepcopy(blocks))
    lists = [
        block["sequence"] for block in blocks if block["blocktype"] == "barcodeList"
    ]
    held = [
        list(decoder.sbcs if isinstance(decoder, SBCDecoder) else decoder.bcs)
        for decoder in decoders
    ]
    assert held == lists, "decoder i must hold list i"
    return [
        "FREE" if isinstance(decoder, SBCDecoder) else "Levenshtein"
        for decoder in decoders
    ]


AGREEING_LAYOUTS = [
    "LCLCLR",
    "LCLC",
    "LCLR",
    "LCL",
    "LL",
    "LRC",
    "LC",
    "LLLC",
    "LLL",
    "L",
    "LR",
    "RL",
    "CL",
    "CRC",
]
DIFFERING_LAYOUTS = ["LRCL", "LRL", "LRLC", "LRCLC"]


@pytest.mark.parametrize("layout", AGREEING_LAYOUTS)
def test_decoders_agree_with_bcwithqc(core, layout):
    blocks = layout_blocks(layout)
    gui = [label for label in core.call("decoder_for_each_list", blocks) if label]
    assert gui == bcwithqc_decoder_labels(blocks)


@pytest.mark.xfail(
    strict=True,
    reason="bcwithqc #decoder-selection: only a read's last list can be FREE",
)
@pytest.mark.parametrize("layout", DIFFERING_LAYOUTS)
def test_decoders_differ_where_bcwithqc_lags_the_rule(core, layout):
    blocks = layout_blocks(layout)
    gui = [label for label in core.call("decoder_for_each_list", blocks) if label]
    assert gui == bcwithqc_decoder_labels(blocks)


# ---------------------------------------------------------------------------
# F2: end to end on error-free reads


def example(name):
    return without_names(read_config(EXAMPLES_DIR / name))


END_TO_END = {
    "sdr_seq_rna": lambda: example("cDNA_paired_end_config.json"),
    "sdr_seq_gdna": lambda: example("gDNA_paired_end_config.json"),
    "pe_mini": lambda: example("pe_mini_config.json"),
    "se_mini": lambda: example("se_mini_config.json"),
    "no_list_on_read_1": lambda: make_config(
        [constant_region("con_r1"), random_block("umi_r1", length=6)],
        [barcode_list("bc_r2"), constant_region("con_r2")],
    ),
    "discarded_last_list": lambda: make_config(
        [
            barcode_list("cb1"),
            constant_region("con1", keep=True),
            barcode_list("cb2", SAMPLE_BARCODES_8, function="discard"),
            random_block("umi", length=6),
        ],
        keep_r1=True,
    ),
}


@pytest.mark.parametrize("name", sorted(END_TO_END))
def test_bcwithqc_output_matches_the_outputs_panel(core, tmp_path, bcwithqc_env, name):
    run_end_to_end(core, END_TO_END[name](), tmp_path, bcwithqc_env)


# ---------------------------------------------------------------------------
# F3: the rules that protect against bcwithqc failures


F3_BASE = make_config(
    [
        barcode_list("cb1"),
        constant_region("con"),
        barcode_list("sb", SAMPLE_BARCODES_8),
        random_block("umi", length=6),
    ],
    keep_r1=True,
)


def f3_config(change):
    config = copy.deepcopy(F3_BASE)
    change(config["barcode_struct_r1"]["blocks"])
    return config


F3_CASES = {
    "E03 empty Levenshtein list": ("E03", lambda b: b[0].update(sequence=[])),
    "E03 empty FREE list": ("E03", lambda b: b[2].update(sequence=[])),
    "E03 empty constant": ("E03", lambda b: b[1].update(sequence=[])),
    "E03 N in a FREE list": (
        "E03",
        lambda b: b[2]["sequence"].__setitem__(0, "TCGCNTTA"),
    ),
    "E05 mixed lengths on bcwithqc's FREE list": (
        "E05",
        lambda b: b[2]["sequence"].__setitem__(0, "TCGCCTTAA"),
    ),
    "E07 zero-length last block": ("E07", lambda b: b[3].update(length=0)),
    "E07 negative last block": ("E07", lambda b: b[3].update(length=-2)),
    "E11 text maxerrors on a list": ("E11", lambda b: b[0].update(maxerrors="1")),
    "E11 fraction maxerrors on a list": ("E11", lambda b: b[0].update(maxerrors=1.5)),
    "E11 text maxerrors on a constant": ("E11", lambda b: b[1].update(maxerrors="1")),
    "E11 text length": ("E11", lambda b: b[3].update(length="8")),
    "E11 unknown type without a length": (
        "E11",
        lambda b: b[1].update(blocktype="linker"),
    ),
}


@pytest.fixture(scope="module")
def f3_reads(tmp_path_factory, bcwithqc_env):
    folder = tmp_path_factory.mktemp("f3_reads")
    return simulate(F3_BASE, folder, bcwithqc_env, nreads=200, seed=3)


def test_f3_base_config_runs(core, tmp_path, bcwithqc_env, f3_reads):
    text = gui_export(core, F3_BASE)["text"]
    _, result = preprocess(text, tmp_path, f3_reads, bcwithqc_env)
    assert result.returncode == 0, result.stderr[-3000:]


@pytest.mark.parametrize("case", sorted(F3_CASES))
def test_rule_guards_against_a_real_bcwithqc_failure(
    core, tmp_path, bcwithqc_env, f3_reads, case
):
    rule, change = F3_CASES[case]
    broken = f3_config(change)
    assert rule in [issue["rule"] for issue in gui_export(core, broken)["issues"]]
    _, result = preprocess(
        json.dumps(broken, indent=4) + "\n", tmp_path, f3_reads, bcwithqc_env
    )
    assert result.returncode != 0, (
        f"bcwithqc accepted a config the GUI rejects with {rule}"
    )


def test_maxerrors_written_as_one_point_zero_breaks_bcwithqc_but_the_gui_writes_one(
    core, tmp_path, bcwithqc_env, f3_reads
):
    """The GUI reads 1.0 as 1 and writes 1, because bcwithqc fails on 1.0."""
    text = gui_export(core, F3_BASE)["text"].replace(
        '"maxerrors": 1,', '"maxerrors": 1.0,', 1
    )
    assert json.loads(text)["barcode_struct_r1"]["blocks"][0]["maxerrors"] == 1.0
    _, result = preprocess(text, tmp_path, f3_reads, bcwithqc_env)
    assert result.returncode != 0


def random_valid_config(rng):
    """A random config the GUI accepts: 1 to 5 blocks on Read 1, short lists, maxerrors up to 2."""
    letters = [rng.choice("LCR") for _ in range(rng.randint(1, 5))]
    blocks = []
    for index, letter in enumerate(letters):
        name = f"b{index}"
        if letter == "L":
            length = rng.randint(6, 9)
            sequences = sorted(
                {
                    "".join(rng.choice("ACGT") for _ in range(length))
                    for _ in range(rng.randint(2, 8))
                }
            )
            blocks.append(barcode_list(name, sequences, maxerrors=rng.randint(0, 2)))
        elif letter == "C":
            variants = list(
                dict.fromkeys(
                    "".join(rng.choice("ACGT") for _ in range(rng.randint(8, 14)))
                    for _ in range(rng.randint(1, 3))
                )
            )
            blocks.append(constant_region(name, variants, maxerrors=rng.randint(0, 2)))
        else:
            blocks.append(
                random_block(
                    name,
                    length=rng.randint(1, 10),
                    function=rng.choice(["UMI", "barcode", "discard"]),
                )
            )
    return make_config(blocks, keep_r1=rng.random() < 0.5)


def test_random_configs_the_gui_accepts_pass_bcwithqcs_builders(core):
    from bcwithqc.misc import build_bc_aligners, build_bc_decoders

    rng = random.Random(20261104)
    configs = [random_valid_config(rng) for _ in range(300)]
    verdicts = core.run(
        "(core, arg) => arg.map((config) => core.validate_config(core.config_from_json(config).model).some((issue) => issue.level === 'error'))",
        configs,
    )
    accepted = [config for config, has_error in zip(configs, verdicts) if not has_error]
    assert len(accepted) > 200
    for config in accepted:
        blocks = config["barcode_struct_r1"]["blocks"]
        build_bc_aligners(copy.deepcopy(blocks), False)
        build_bc_decoders(copy.deepcopy(blocks))


# ---------------------------------------------------------------------------
# F4: temporary warnings. Each case shows today's failure; once it passes,
# delete the warning from TEMPORARY_WARNINGS in the page, and this case.


@pytest.mark.xfail(
    strict=True, reason="bcwithqc #one-block-crash: IndexError on a one-block read"
)
def test_w11_one_block_read(core, tmp_path, bcwithqc_env):
    run_end_to_end(
        core,
        make_config([barcode_list("bc", SAMPLE_BARCODES_8)], keep_r1=True),
        tmp_path,
        bcwithqc_env,
    )


@pytest.mark.xfail(
    strict=True,
    reason="bcwithqc #r2-only-crash: Read 1 without blocks crashes the aligner",
)
def test_w12_barcodes_on_read_2_only(core, tmp_path, bcwithqc_env):
    """Counts as fixed only when everything matches, the exact CB string included (decision D13)."""
    config = make_config(
        [],
        [
            barcode_list("bc_r2"),
            constant_region("con_r2"),
            random_block("umi_r2", length=6),
        ],
    )
    run_end_to_end(core, config, tmp_path, bcwithqc_env)


@pytest.mark.xfail(
    strict=True,
    reason="bcwithqc #barcode-file-guess: without a constant region it picks file 2",
)
def test_w08_barcode_file_guess(core, tmp_path, bcwithqc_env):
    config = make_config(
        [barcode_list("bc", SAMPLE_BARCODES_8), random_block("umi", length=6)], []
    )
    exported = gui_export(core, config)
    # Long Read 2 reads make bcwithqc's guess pick the wrong file.
    reads = simulate(
        json.loads(exported["text"]),
        tmp_path,
        bcwithqc_env,
        nreads=400,
        extra=("--random-tail-length=60",),
    )
    output, result = preprocess(exported["text"], tmp_path, reads, bcwithqc_env)
    assert result.returncode == 0, result.stderr[-2000:]
    # bcwithqc names the tags file after the file it takes the barcodes from.
    first_fastq = fastq_files(reads)[0]
    prefix = first_fastq.name.removesuffix(".txt.gz")
    tags_files = [path.name for path in output.glob("rec_names_and_tags1_*.txt")]
    assert tags_files == [f"rec_names_and_tags1_{prefix}.txt"], (
        "barcodes were read from the wrong file"
    )
    first_file = {name.split("_")[0]: name for name, _ in read_fastq(first_fastq)}
    blocks = json.loads(exported["text"])["barcode_struct_r1"]["blocks"]
    lines = (output / tags_files[0]).read_text().splitlines()
    assert len(lines) == 400
    for line in lines:
        name, *fields = line.split("\t")
        tags = dict(field.split(":", 2)[::2] for field in fields)
        values, _ = truth(first_file[name.split("_")[0]], blocks)
        assert tags["CB"] == values["bc"]


@pytest.mark.xfail(
    strict=True, reason="bcwithqc #orientation-cap: NameError on read pair 20,011"
)
def test_w05_orientation_cap(core, tmp_path, bcwithqc_env):
    config = make_config(
        [barcode_list("bc"), constant_region("con"), random_block("umi", length=6)],
        [],
        orientation=True,
    )
    exported = gui_export(core, config)
    plain = json.loads(exported["text"])
    plain["unknown_read_orientation"] = False
    reads = simulate(
        plain, tmp_path, bcwithqc_env, nreads=20_011, extra=("--random-tail-length=30",)
    )
    output, result = preprocess(exported["text"], tmp_path, reads, bcwithqc_env)
    assert result.returncode == 0, result.stderr[-2000:]
    (tags_file,) = output.glob("rec_names_and_tags1_*.txt")
    assert len(tags_file.read_text().splitlines()) == 20_011


def test_sample_constant_is_the_examples_constant():
    """F3_BASE and the F4 cases use the pe_mini constant, so reads match bcwithqc's own examples."""
    blocks = read_config(EXAMPLES_DIR / "pe_mini_config.json")["barcode_struct_r1"][
        "blocks"
    ]
    assert blocks[1]["sequence"] == [SAMPLE_CONSTANT]
    assert blocks[0]["sequence"] == SAMPLE_BARCODES


@pytest.mark.xfail(
    strict=True,
    reason="bcwithqc #decoder-picks-farthest: DistanceThresh negates pywfa's positive distance",
)
def test_w14_levenshtein_picks_the_farthest_barcode():
    from bcwithqc.bc_decoders import BCDecoder

    nearest, farther = "ACGTACGTA", "ACGTTGCTA"
    read = "ACGTTCGTA"  # 1 edit from the first, 2 from the second
    assert BCDecoder([nearest, farther], 2).decode(read) == nearest


# ---------------------------------------------------------------------------
# Outputs that bcwithqc writes differently from the panel today
# (TEMPORARY_OUTPUT_NOTES in the page).


@pytest.mark.xfail(
    strict=True,
    reason="bcwithqc #no-umi-outputs: UR is written empty without a UMI block",
)
def test_no_umi_block_writes_no_ur_tag(core, tmp_path, bcwithqc_env):
    config = make_config(
        [
            barcode_list("cb1"),
            constant_region("con1"),
            barcode_list("cb2", SAMPLE_BARCODES_8),
        ],
        keep_r1=True,
    )
    exported = gui_export(core, config)
    assert exported["outputs"]["UR"]["written"] is False
    reads = simulate(json.loads(exported["text"]), tmp_path, bcwithqc_env, nreads=200)
    output, result = preprocess(exported["text"], tmp_path, reads, bcwithqc_env)
    assert result.returncode == 0, result.stderr[-2000:]
    assert count(tmp_path, output, bcwithqc_env).returncode == 0
    records = bam_tags(output / "with_bc_umi.sorted.bam")
    assert records
    assert not [name for name, tags in records if "UR" in tags]


@pytest.mark.xfail(
    strict=True,
    reason="bcwithqc #empty-alignment-files: a read that keeps nothing still gets a file",
)
def test_read_that_keeps_nothing_gets_no_alignment_file(core, tmp_path, bcwithqc_env):
    config = make_config(
        [barcode_list("cb1"), constant_region("con1"), random_block("umi", length=6)],
        [],
    )
    exported = gui_export(core, config)
    assert exported["outputs"]["alignment_files"][0]["kind"] == "none"
    reads = simulate(json.loads(exported["text"]), tmp_path, bcwithqc_env, nreads=100)
    output, result = preprocess(exported["text"], tmp_path, reads, bcwithqc_env)
    assert result.returncode == 0, result.stderr[-2000:]
    first_prefix = fastq_files(reads)[0].name.removesuffix(".txt.gz")
    assert not (output / f"sans_bc_{first_prefix}.fq").exists()


# ---------------------------------------------------------------------------
# D12: mistakes bcwithqc accepts silently today. Each should be rejected once
# bcwithqc validates configs on load (notes #config-validation).

SILENTLY_ACCEPTED = {
    "unknown block type with a length": lambda b: b[3].update(blocktype="umiBlock"),
    "unknown block function": lambda b: b[3].update(blockfunction="umi"),
    "text maxerrors on a random block": lambda b: b[3].update(maxerrors="1"),
    "zero-length random block that is not last": lambda b: b.insert(
        1, random_block("gap", length=0)
    ),
}


@pytest.mark.xfail(
    strict=True, reason="bcwithqc #config-validation: configs are not checked on load"
)
@pytest.mark.parametrize("case", sorted(SILENTLY_ACCEPTED))
def test_bcwithqc_rejects_what_the_gui_rejects(
    core, tmp_path, bcwithqc_env, f3_reads, case
):
    broken = f3_config(SILENTLY_ACCEPTED[case])
    assert [
        issue
        for issue in gui_export(core, broken)["issues"]
        if issue["level"] == "error"
    ]
    _, result = preprocess(
        json.dumps(broken, indent=4) + "\n", tmp_path, f3_reads, bcwithqc_env
    )
    assert result.returncode != 0


# ---------------------------------------------------------------------------
# D14: further cases


def all_layouts(longest):
    """Every layout of 1 to `longest` blocks, as strings of L, C and R."""
    return [
        "".join(letters)
        for length in range(1, longest + 1)
        for letters in itertools.product("LCR", repeat=length)
    ]


def bcwithqc_lags_the_rule(layout):
    """Derived from bcwithqc's code: a list other than the read's last one is followed by a random block."""
    last_list = layout.rfind("L")
    return any(
        letter == "L" and index != last_list and layout[index + 1 : index + 2] == "R"
        for index, letter in enumerate(layout)
    )


EVERY_LAYOUT = [
    pytest.param(
        layout,
        marks=pytest.mark.xfail(strict=True, reason="bcwithqc #decoder-selection"),
    )
    if bcwithqc_lags_the_rule(layout)
    else layout
    for layout in all_layouts(5)
]


def test_every_layout_list_has_363_entries():
    assert len(all_layouts(5)) == 3 + 9 + 27 + 81 + 243


@pytest.mark.parametrize("layout", EVERY_LAYOUT)
def test_decoders_agree_for_every_layout_up_to_five_blocks(core, layout):
    blocks = layout_blocks(layout)
    gui = [label for label in core.call("decoder_for_each_list", blocks) if label]
    assert gui == bcwithqc_decoder_labels(blocks)


@pytest.mark.xfail(
    strict=True,
    reason="bcwithqc #discard-slip: a later list reports the discarded list's barcode",
)
def test_w01_discarded_list_shifts_later_barcodes(core, tmp_path, bcwithqc_env):
    config = make_config(
        [
            barcode_list("cb1", function="discard"),
            constant_region("con1", keep=True),
            barcode_list("cb2", SAMPLE_BARCODES_8),
            random_block("umi", length=6),
        ],
        keep_r1=True,
    )
    run_end_to_end(core, config, tmp_path, bcwithqc_env)


def test_mate_swapped_pairs_are_decoded_with_orientation_on(
    core, tmp_path, bcwithqc_env
):
    """Below the 20,010-pair cap, orientation on finds the barcodes in either file."""
    config = example("cDNA_paired_end_config.json")
    config["unknown_read_orientation"] = True
    exported = gui_export(core, config)
    plain = json.loads(exported["text"])
    plain["unknown_read_orientation"] = False
    reads = simulate(plain, tmp_path, bcwithqc_env, nreads=600)
    first, second = (read_fastq(path) for path in fastq_files(reads))
    swapped = tmp_path / "swapped"
    swapped.mkdir()
    out_first, out_second = [], []
    for index, (record_1, record_2) in enumerate(zip(first, second)):
        pair = (record_2, record_1) if index % 3 == 0 else (record_1, record_2)
        out_first.append(pair[0])
        out_second.append(pair[1])
    import gzip

    for name, records in (
        ("swapped_1.txt.gz", out_first),
        ("swapped_2.txt.gz", out_second),
    ):
        with gzip.open(swapped / name, "wt") as handle:
            for read_name, sequence in records:
                handle.write(f"@{read_name}\n{sequence}\n+\n{'E' * len(sequence)}\n")
    output, result = preprocess(exported["text"], tmp_path, swapped, bcwithqc_env)
    assert result.returncode == 0, result.stderr[-2000:]
    (tags_file,) = output.glob("rec_names_and_tags1_*.txt")
    blocks = plain["barcode_struct_r1"]["blocks"]
    truth_by_name = {name: truth(name, blocks)[0] for name, _ in first}
    lines = tags_file.read_text().splitlines()
    assert len(lines) == len(first)
    for line in lines:
        name, *fields = line.split("\t")
        tags = dict(field.split(":", 2)[::2] for field in fields)
        values = truth_by_name[name]
        assert tags["CB"] == ".".join(values[key] for key in ("cb1", "cb2", "sb"))


def test_slash_in_a_block_name_makes_count_fail(core, tmp_path, bcwithqc_env, f3_reads):
    """E08 blocks "/" in names (decision D7) because bcwithqc fails on it."""
    broken = f3_config(lambda b: b[0].update(blockname="cb/1"))
    assert "E08" in [issue["rule"] for issue in gui_export(core, broken)["issues"]]
    output, result = preprocess(
        json.dumps(broken, indent=4) + "\n", tmp_path, f3_reads, bcwithqc_env
    )
    if result.returncode == 0:
        result = count(tmp_path, output, bcwithqc_env)
    assert result.returncode != 0
