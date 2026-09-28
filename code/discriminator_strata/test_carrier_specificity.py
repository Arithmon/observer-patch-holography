from __future__ import annotations

import pytest

import carrier_specificity as cs


def test_every_member_builds_with_declared_counts() -> None:
    expected = {
        "icosahedron": (12, 30),
        "hemi_icosahedron": (6, 15),
        "cuboctahedron": (12, 24),
        "truncated_tetrahedron": (12, 18),
        "hexagonal_prism": (12, 18),
        "octahedron": (6, 12),
        "cube": (8, 12),
        "dodecahedron": (20, 30),
        "tetrahedron": (4, 6),
        "petersen": (10, 15),
        "complete_k12": (12, 66),
        "cycle_c12": (12, 12),
        "complete_bipartite_k66": (12, 36),
    }
    for name, (vertices, edges) in expected.items():
        row = cs.probe_member(name)
        assert (row["vertices"], row["edges"]) == (vertices, edges), name


def test_icosahedron_is_the_unique_full_hit() -> None:
    receipt = cs.build_receipt()
    assert receipt["full_hit_members"] == ["icosahedron"]
    assert receipt["unique_full_hit"] is True
    assert receipt["ensemble_size"] == 13


def test_cuboctahedron_shares_multiplicities_and_fails_elsewhere() -> None:
    row = cs.probe_member("cuboctahedron")
    assert row["band_multiplicities"] == [1, 3, 3, 5]
    assert row["probes"]["multiplicities_1_3_3_5"] is True
    assert row["probes"]["degree_five_regular"] is False
    assert row["probes"]["galois_band_pair"] is False
    assert row["probes"]["response_involution"] is False
    assert row["full_hit"] is False


def test_per_probe_counts_are_stable() -> None:
    receipt = cs.build_receipt()
    assert receipt["per_probe_pass_counts"] == {
        "twelve_ports": 7,
        "degree_five_regular": 2,
        "four_bands": 3,
        "galois_band_pair": 3,
        "multiplicities_1_3_3_5": 2,
        "unique_antipode_involution": 7,
        "twenty_triangles": 2,
        "response_involution": 1,
    }


def test_galois_pair_members() -> None:
    for name, expected in (
        ("icosahedron", True),
        ("dodecahedron", True),
        ("cycle_c12", True),
        ("petersen", False),
        ("cube", False),
    ):
        row = cs.probe_member(name)
        assert row["probes"]["galois_band_pair"] is expected, name


def test_response_law_family_counts_and_declared_recovery() -> None:
    laws = cs.response_law_family()
    assert laws["law_count"] == 16
    assert laws["rational_law_count"] == 8
    assert laws["declared_coefficients_recovered"] is True
    declared = [row for row in laws["laws"] if row["declared_law"]]
    assert len(declared) == 1
    assert declared[0]["rational_coefficients"] is True
    assert declared[0]["coefficients_low_to_high"] == [
        "1+0*sqrt5",
        "-1/2+0*sqrt5",
        "-2/5+0*sqrt5",
        "1/10+0*sqrt5",
    ]
    for row in laws["laws"]:
        symmetric = row["sign_pattern"][1] == row["sign_pattern"][2]
        assert row["rational_coefficients"] is symmetric


def test_hemi_icosahedron_is_the_antipodal_quotient() -> None:
    edges = cs.hemi_icosahedron()
    names = sorted({name for edge in edges for name in edge})
    assert len(names) == 6
    assert sorted(edges) == [
        (a, b) for i, a in enumerate(names) for b in names[i + 1 :]
    ]
    row = cs.probe_member("hemi_icosahedron")
    assert row["band_multiplicities"] == [1, 5]
    hits = {probe for probe, value in row["probes"].items() if value}
    assert hits == {"degree_five_regular", "twenty_triangles"}
    assert (
        row["cell_typing"]["response_involution"]
        == "failed_precondition_non_twelve_port"
    )


def test_band_projectors_are_exact_and_complete() -> None:
    matrix = cs.adjacency_matrix(cs.icosahedron())
    projectors = cs.band_projectors(matrix)
    size = len(matrix)
    total = [
        [sum(p[i][j] for _, p in projectors) for j in range(size)]
        for i in range(size)
    ]
    assert total == [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    for label, p in projectors:
        square = [
            [sum(p[i][k] * p[k][j] for k in range(size)) for j in range(size)]
            for i in range(size)
        ]
        assert square == p, label


def test_galois_pair_is_the_antipode_odd_sector() -> None:
    stratum = cs.orientation_stratum()
    assert stratum["galois_pair_parity"] == "odd"
    assert stratum["even_sector_bands"] == {"-1": 5, "5": 1}
    assert stratum["odd_sector_bands"] == {"galois_pair": 6}
    assert stratum["quotient_bands"] == stratum["even_sector_bands"]
    assert stratum["quotient_bands_equal_even_sector"] is True
    assert stratum["charge_mode_band"] == "5"
    assert stratum["charge_mode_parity"] == "even"


def test_replayed_dimension_argument_admits_nothing_on_the_quotient() -> None:
    twelve = cs.replay_dimension_argument(12)
    assert twelve["candidates"] == [
        "su(2)+su(2)+su(2)+su(2)",
        "u(1)+su(2)+su(3)",
    ]
    assert twelve["admissible"] == ["u(1)+su(2)+su(3)"]
    six = cs.replay_dimension_argument(6)
    assert six["candidates"] == ["su(2)+su(2)"]
    assert six["admissible"] == []
    with pytest.raises(cs.SpecificityError, match="twelve ports"):
        cs.replay_dimension_argument(13)


def test_antipode_parity_table_marks_every_galois_pair_odd() -> None:
    table = cs.antipode_parity_table()
    assert sorted(table) == [
        "cube",
        "cuboctahedron",
        "cycle_c12",
        "dodecahedron",
        "hexagonal_prism",
        "icosahedron",
        "octahedron",
    ]
    galois = {
        name: [row["parity"] for row in rows if row["band"] == "galois_pair"]
        for name, rows in table.items()
    }
    for name in ("icosahedron", "dodecahedron", "cycle_c12"):
        assert galois[name] == ["odd"], name
    for name in ("cube", "cuboctahedron", "hexagonal_prism", "octahedron"):
        assert galois[name] == [], name


def test_committed_receipt_is_byte_exact() -> None:
    committed = cs.OUTPUT_PATH.read_bytes()
    assert committed == cs.canonical_json_bytes(cs.build_receipt())


def test_duplicate_edge_fails_closed() -> None:
    with pytest.raises(cs.SpecificityError, match="duplicate edge"):
        cs.adjacency_matrix([("a", "b"), ("b", "a")])


def test_minimal_polynomial_on_known_graph() -> None:
    matrix = cs.adjacency_matrix(cs.complete_k12())
    minimal = cs.minimal_polynomial(matrix)
    assert len(minimal) - 1 == 2
    structure = cs.spectral_structure(matrix)
    assert structure["rational_bands"] == ["-1", "11"]
