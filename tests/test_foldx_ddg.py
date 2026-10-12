"""Tests for the FoldX helpers that do not need the FoldX binary."""
import pytest

from src.structure.foldx_ddg import parse_dif, shard_chunks

DIF = (
    "FoldX 5.1 (c)\n"
    "PDB file analysed: X_Repair.pdb\n"
    "Output type: BuildModel\n"
    "Pdb\ttotal energy\tBackbone Hbond\n"
    "X_Repair_1.pdb\t1.25\t0.1\n"
    "X_Repair_2.pdb\t-0.50\t0.0\n"
)


def test_parse_dif_returns_total_energy_in_order(tmp_path):
    path = tmp_path / "Dif_X_Repair.fxout"
    path.write_text(DIF)
    assert parse_dif(path, 2, runs=1) == [1.25, -0.5]
    assert parse_dif(path, 1, runs=2) == [0.375]


def test_parse_dif_rejects_wrong_row_count(tmp_path):
    path = tmp_path / "Dif_X_Repair.fxout"
    path.write_text(DIF)
    with pytest.raises(ValueError):
        parse_dif(path, 3, runs=1)


def test_shards_cover_every_chunk_once_and_follow_weights():
    chunks = [{"chunk": f"c{i}", "n_variants": 100, "n_residues": 300} for i in range(9)]
    shards = [shard_chunks(chunks, [2, 1], k) for k in range(2)]
    assert shards[0] | shards[1] == {c["chunk"] for c in chunks}
    assert not shards[0] & shards[1]
    assert len(shards[0]) == 6 and len(shards[1]) == 3
