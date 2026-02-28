from hypnotoad.agent import tools


def test_validate_unknown_key():
    # Using an invalid key should lead to an issue and suggestion
    assert tools.validate_settings({"nxcore": 10}) == {
        "valid": False,
        "issues": {
            "nxcore": {
                "issue_type": "unknown_key",
                "message": "'nxcore' is not a recognised Hypnotoad option",
                "suggestions": ["nx_core"],
            }
        },
    }


def test_validate_wrong_type():
    # Using the wrong type
    assert tools.validate_settings({"nx_core": 3.4}) == {
        "valid": False,
        "issues": {
            "nx_core": {
                "issue_type": "wrong_type",
                "message": "'nx_core' expects int, got float",
                "expected": "int",
                "got": "float",
            }
        },
    }


def test_validate_invalid_value():
    assert tools.validate_settings({"curvature_type": "nonsense"}) == {
        "valid": False,
        "issues": {
            "curvature_type": {
                "issue_type": "invalid_value",
                "message": "'curvature_type' = 'nonsense' not in allowed values: ('curl(b/B)', 'curl(b/B) with x-y derivatives', 'bxkappa')",
                "expected": "('curl(b/B)', 'curl(b/B) with x-y derivatives', 'bxkappa')",
                "got": "nonsense",
            }
        },
    }


def test_validate_out_of_bounds():
    assert tools.validate_settings({"refine_width": -1.0}) == {
        "valid": False,
        "issues": {
            "refine_width": {
                "issue_type": "out_of_bounds",
                "message": "'refine_width' = -1.0 failed validation: Failed check: def is_positive(x):\n    try:\n        return x > 0\n    except TypeError:\n        return False",
                "got": "-1.0",
            }
        },
    }
