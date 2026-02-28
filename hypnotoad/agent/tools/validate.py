"""
Validate a dict of Hypnotoad settings, returning a dict
describing issues.
"""

import difflib
import inspect
from dataclasses import dataclass
from enum import Enum
from typing import Any

from ...cases import tokamak
from ...core.mesh import BoutMesh

POSSIBLE_OPTIONS = (
    tokamak.TokamakEquilibrium.user_options_factory.defaults
    | tokamak.TokamakEquilibrium.nonorthogonal_options_factory.defaults
    | BoutMesh.user_options_factory.defaults
)

POSSIBLE_KEYS = [opt for opt in POSSIBLE_OPTIONS]


class IssueType(str, Enum):
    UNKNOWN_KEY = "unknown_key"
    WRONG_TYPE = "wrong_type"
    OUT_OF_BOUNDS = "out_of_bounds"
    INVALID_VALUE = "invalid_value"
    MISSING = "missing"


@dataclass
class SettingIssue:
    issue_type: IssueType
    message: str
    suggestions: list[str] = None
    expected: Any = None
    got: Any = None

    def to_dict(self) -> dict:
        d = {"issue_type": self.issue_type.value, "message": self.message}
        if self.suggestions:
            d["suggestions"] = self.suggestions
        if self.expected is not None:
            d["expected"] = str(self.expected)
        if self.got is not None:
            d["got"] = str(self.got)
        return d


def _describe_type(value_type) -> str:
    """Human-readable type description from value_type field."""
    if value_type is None:
        return "any"
    if isinstance(value_type, (list, tuple)):
        names = [t.__name__ if t is not None else "None" for t in value_type]
        return " or ".join(names)
    return value_type.__name__


def _check_value_type(value: Any, value_type) -> bool:
    """Check value against value_type, which may be a type, list of types, or None."""
    if value_type is None:
        return True  # Any type accepted
    types = value_type if isinstance(value_type, (list, tuple)) else [value_type]
    # Treat None in the type list as allowing Python None
    allowed_types = tuple(t for t in types if t is not None)
    none_allowed = any(t is None for t in types)
    if value is None:
        return none_allowed
    # Be lenient: int is acceptable where float is expected
    if isinstance(value, bool):
        # bool is a subclass of int -- only allow if bool is explicitly in types
        return bool in types or (not allowed_types)
    if isinstance(value, int) and float in types:
        return True
    return isinstance(value, allowed_types)


def _run_checks(value: Any, meta: WithMeta) -> list[str]:
    """
    Run check_all and check_any against value.
    Returns list of failure messages (empty = all passed).
    """
    failures = []

    # check_all: every check must pass
    check_all = meta.check_all if hasattr(meta, "check_all") else None
    if check_all is not None:
        checks = check_all if isinstance(check_all, (list, tuple)) else [check_all]
        for check in checks:
            try:
                if not check(value):
                    # Try to get a useful description of the check
                    src = (
                        inspect.getsource(check).strip()
                        if hasattr(check, "__code__")
                        else str(check)
                    )
                    failures.append(f"Failed check: {src}")
            except Exception as e:
                failures.append(f"Check raised exception: {e}")

    # check_any: at least one check must pass
    check_any = meta.check_any if hasattr(meta, "check_any") else None
    if check_any is not None:
        checks = check_any if isinstance(check_any, (list, tuple)) else [check_any]
        passed = False
        for check in checks:
            try:
                if check(value):
                    passed = True
                    break
            except Exception:
                pass
        if not passed:
            srcs = []
            for check in checks:
                try:
                    srcs.append(inspect.getsource(check).strip())
                except Exception:
                    srcs.append(str(check))
            failures.append(f"Must satisfy at least one of: {'; '.join(srcs)}")

    return failures


def validate_settings(settings: dict = {}) -> dict:
    """Check settings for common issues before running.

    Returns a dictionary with a boolean flag 'valid'
    and a dict of 'issues' indexed by keys that do not
    pass validation.
    """
    issues = {}

    for key, value in settings.items():
        # Unknown keys
        if key not in POSSIBLE_KEYS:
            suggestions = difflib.get_close_matches(key, POSSIBLE_KEYS, n=3, cutoff=0.6)
            issues[key] = SettingIssue(
                issue_type=IssueType.UNKNOWN_KEY,
                message=f"'{key}' is not a recognised Hypnotoad option",
                suggestions=suggestions or None,
            ).to_dict()
            continue

        # Validate using WithMeta
        meta = POSSIBLE_OPTIONS[key]
        value_type = getattr(meta, "value_type", None)
        allowed = getattr(meta, "allowed", None)

        if not _check_value_type(value, value_type):
            issues[key] = SettingIssue(
                issue_type=IssueType.WRONG_TYPE,
                message=(
                    f"'{key}' expects {_describe_type(value_type)}, "
                    f"got {type(value).__name__}"
                ),
                expected=_describe_type(value_type),
                got=type(value).__name__,
            ).to_dict()
            continue  # no point checking further with wrong type

        # Allowed values check
        if allowed is not None and value not in allowed:
            suggestions = difflib.get_close_matches(
                str(value), [str(v) for v in allowed], n=3, cutoff=0.5
            )
            issues[key] = SettingIssue(
                issue_type=IssueType.INVALID_VALUE,
                message=f"'{key}' = {value!r} not in allowed values: {allowed}",
                expected=allowed,
                got=value,
                suggestions=suggestions or None,
            ).to_dict()
            continue

        # check_all / check_any
        check_failures = _run_checks(value, meta)
        if check_failures:
            issues[key] = SettingIssue(
                issue_type=IssueType.OUT_OF_BOUNDS,
                message=f"'{key}' = {value!r} failed validation: {'; '.join(check_failures)}",
                got=value,
            ).to_dict()

    return {"valid": len(issues) == 0, "issues": issues}
