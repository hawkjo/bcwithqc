# Sources of the GUI test fixtures

The config builder's tests compare its output with configs written by other
tools. This file records where each fixture came from, by commit.

## Copied files

| File | Source | Notes |
| --- | --- | --- |
| `legacy_master_cDNA.json` | bcwithqc `examples/cDNA.json` at `6b36b7b1ff75d90156b0e4d88dda76c0a56d4299` | Config format from before `single_end_reads` and `keepblock`. Byte-identical to SDRranger `examples/cDNA.json` at `548683d7ac2467f2d65d3b44d4a1756a4d2d47d5`. |
| `legacy_master_gDNA.json.gz` | bcwithqc `examples/gDNA.json.gz` at `6b36b7b1ff75d90156b0e4d88dda76c0a56d4299` | Same format, gzipped. |
| `legacy_sdrranger_gDNA.json.gz` | SDRranger (https://github.com/hawkjo/SDRranger) `examples/gDNA.json.gz` at `548683d7ac2467f2d65d3b44d4a1756a4d2d47d5` | Same format; `keep_nonbarcode` is `false` on Read 1, unlike the bcwithqc copy. |

## Files read in place

The current-format configs in `examples/*.json` are read from the checkout
under test rather than copied, so the tests follow any change to them. At
`5d08d15d0fdb803fd89d87d48d4f0aa68241f92a` (the merge of PR #8) they are the
same as on `reopen-lukas-branch`.

The example config in `doc/config_format.md` is also read in place: the tests
extract the JSON block under the "Example Config.json" heading.

## Generated files

The old barcodinggui `.json.gz` used by the migration test is generated at
test time by `write_old_gui_config` in `gui_support.py`. It follows
barcodinggui 2b473d1 (`mainwindow.MainWindow._save`,
`readprops.ReadProperties.to_dict`, `sequenceprops.SequenceProps.to_dict`):
`json.dump(config, indent=4)` into `gzip.open(..., "wt", compresslevel=9)`,
default block names `sequence 0`, `sequence 1`, ..., a `name` list on every
constant region, and no `keepblock` or `single_end_reads` keys.
