def test_compute_formation_value_for_known_compound(ops_context):
    value = ops_context.thermo_ops.compute_formation_value(
        ops_context.water["cid"],
        20.0,
        {
            101: 2.0,
            102: 8.0,
            103: 1.0,
        },
    )

    assert value == 14.0


def test_compute_formation_value_returns_none_for_unknown_or_missing_elements(ops_context):
    assert ops_context.thermo_ops.compute_formation_value(9999, 10.0, {101: 1.0}) is None
    assert (
        ops_context.thermo_ops.compute_formation_value(
            ops_context.water["cid"],
            10.0,
            {101: 1.0},
        )
        is None
    )
