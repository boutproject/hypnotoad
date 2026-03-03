from dataclasses import dataclass
from typing import Optional
import inspect


@dataclass
class Chunk:
    text: str
    section: str
    source: str
    chunk_type: str  # "rst", "function", "usage", "option"
    name: Optional[str] = None  # function/option name if applicable
    start_line: Optional[int] = None
    url: Optional[str] = None

    def to_dict(self) -> dict:
        d = {
            "section": self.section,
            "source": self.source,
            "chunk_type": self.chunk_type,
            "text": self.text,
        }
        if self.name:
            d["name"] = self.name
        if self.start_line:
            d["start_line"] = self.start_line
        if self.url:
            d["url"] = self.url
        return d


def _describe_type(value_type) -> str:
    """Human-readable type description from value_type field."""
    if value_type is None:
        return "any"
    if isinstance(value_type, (list, tuple)):
        names = [t.__name__ if t is not None else "None" for t in value_type]
        return " or ".join(names)
    return value_type.__name__


def extract_option_chunks(known_options) -> list[Chunk]:
    """
    Build RAG chunks from the live OptionsFactory.
    Each option becomes one chunk with its full metadata.
    Called once at index-build time.
    """
    chunks = []

    for name, meta in known_options.items():
        default = meta.value  # first arg to WithMeta
        doc = getattr(meta, "doc", "No description available")
        value_type = getattr(meta, "value_type", None)
        allowed = getattr(meta, "allowed", None)
        check_all = getattr(meta, "check_all", None)
        check_any = getattr(meta, "check_any", None)

        # Build a human-readable description for embedding
        lines = [
            f"Option: {name}",
            f"Default: {default}",
            f"Description: {doc}",
        ]
        if value_type is not None:
            lines.append(f"Type: {_describe_type(value_type)}")
        if allowed is not None:
            lines.append(f"Allowed values: {allowed}")
        if check_all is not None:
            checks = check_all if isinstance(check_all, (list, tuple)) else [check_all]
            for c in checks:
                try:
                    lines.append(f"Constraint (all): {inspect.getsource(c).strip()}")
                except Exception:
                    lines.append(f"Constraint (all): {c}")
        if check_any is not None:
            checks = check_any if isinstance(check_any, (list, tuple)) else [check_any]
            for c in checks:
                try:
                    lines.append(f"Constraint (any): {inspect.getsource(c).strip()}")
                except Exception:
                    lines.append(f"Constraint (any): {c}")

        chunks.append(
            Chunk(
                text="\n".join(lines),
                section=f"option: {name}",
                source="hypnotoad.options_factory (live)",
                chunk_type="option",
                name=name,
            )
        )

    return chunks


class ChunkDatabase:
    """
    Index chunks and retrieve based on queries
    """

    def __init__(self, chunks: list[Chunk]):
        from rank_bm25 import BM25Okapi

        self.chunks = chunks
        corpus = [chunk.text.split() for chunk in chunks]
        self.bm25 = BM25Okapi(corpus)

    def retrieve(self, query: str, k: int = 4) -> list[dict]:
        scores = self.bm25.get_scores(query.split())
        top_k = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)[:k]
        return [self.chunks[i].to_dict() for i in top_k]
