from dataclasses import dataclass, asdict
from typing import Optional
import inspect
import numpy as np
from pathlib import Path
import faiss
from rank_bm25 import BM25Okapi
import json


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
    Index chunks and retrieve based on queries.
    This database uses BM25 to rank based on keywords.
    """

    def __init__(self, chunks: list[Chunk]):
        self.chunks = chunks
        corpus = [chunk.text.split() for chunk in chunks]
        self.bm25 = BM25Okapi(corpus)

    def retrieve(self, query: str, k: int = 4) -> list[dict]:
        scores = self.bm25.get_scores(query.split())
        top_k = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)[:k]
        return [self.chunks[i].to_dict() for i in top_k]


class ChunkFaissDatabase:
    """
    Persistent vector database for `Chunk` objects using OpenAI embeddings
    and a FAISS similarity index.

    Each chunk's embedding is stored in a FAISS index (IndexFlatIP) with
    L2-normalised vectors so that inner product corresponds to cosine similarity.

    IMPORTANT:
    The implementation assumes that the order of vectors in the FAISS index
    matches the order of `self.chunks`. The i-th embedding in the index
    corresponds to `self.chunks[i]`. The database is therefore append-only
    and does not support deletion or reordering of chunks without rebuilding
    the index.

    Data is persisted to disk as:
      - faiss.index      : serialized FAISS index
      - manifest.json    : metadata (embedding model, dimension, etc.)
      - chunks.jsonl     : one JSON-encoded Chunk per line

    """

    def __init__(
        self,
        client,
        model: Optional[str] = None,
        restore: Optional[Path | str] = None,
    ):
        """
        Parameters
        ----------
        client : openai.OpenAI
            OpenAI client used to generate embeddings.

        model : str, optional
            Name of embedding model used to generate vectors.
            Required when creating a new database.

        restore : str or Path, optional
            Directory containing previously saved database files.
            If provided, the index and chunks are restored from disk.
        """
        self.client = client
        self.model = model
        self.chunks = []
        self.index = None
        if restore is None:
            # Start a new index
            if model is None:
                raise ValueError("Specify a model for new database")
            # Create index in add_chunks when embedding size is known
        else:
            # Restore database from file
            self.load(restore)
            if (model is not None) and (model != self.model):
                raise ValueError(
                    f"Model '{model}' not equal to '{self.model}' in restored database {restore}"
                )

    def add_chunks(self, chunks: list[Chunk]):
        """
        Add new chunks to the database.

        For each chunk:
        1. Compute its embedding using the configured embedding model.
        2. L2-normalise the embedding (for cosine similarity search).
        3. Append the embedding to the FAISS index.
        4. Append the chunk to `self.chunks`.

        The order of addition is preserved, so the FAISS vector at position i
        corresponds to `self.chunks[i]`.

        Notes
        -----
        - This database is append-only. Removing or reordering chunks will
        break alignment between the FAISS index and `self.chunks`.
        - All embeddings must have the same dimension as the existing index.
        - If the index has not yet been created, it will be initialised
        using the embedding dimension of the first batch.

        Parameters
        ----------
        chunks : list[Chunk]
            Chunks to embed and add to the database.

        Raises
        ------
        ValueError
            If the embedding dimension does not match the existing index.
        """
        if not chunks:
            return

        texts = [chunk.text for chunk in chunks]
        embeddings = self._calculate_embeddings(texts)
        X = np.array(embeddings, dtype="float32")
        faiss.normalize_L2(X)  # for cosine similarity

        if self.index is None:
            self.index = faiss.IndexFlatIP(X.shape[1])
        elif self.index.d != X.shape[1]:
            raise ValueError(f"Embedding dim {X.shape[1]} != index dim {self.index.d}")

        # Add to index and chunks list in the same order
        # so that indices remain synchronised.
        self.index.add(X)
        self.chunks.extend(chunks)

    def retrieve(self, query: str, k: int = 4):
        """
        Retrieve the top-k most similar chunks for a query string.

        The query is embedded using the configured embedding model,
        L2-normalised, and searched against the FAISS index using
        inner product similarity (equivalent to cosine similarity).

        Parameters
        ----------
        query : str
            Natural-language search query.
        k : int, default=4
            Number of top results to return. If k exceeds the number
            of indexed chunks, it will be clamped.

        Returns
        -------
        chunks : list[Chunk]
            Retrieved chunks in descending similarity order.
        scores : list[float]
            Corresponding similarity scores (cosine similarity).

        Notes
        -----
        - Returns empty lists if the index is empty.
        - Scores are inner products of L2-normalised vectors
        (i.e., cosine similarity in [-1, 1]).
        """
        if self.index is None or len(self.chunks) == 0:
            # Nothing added to database
            return [], []
        k = min(k, len(self.chunks))

        query_embedding = self._calculate_embeddings([query])
        q = np.array(query_embedding, dtype="float32")
        faiss.normalize_L2(q)

        Dists, Inds = self.index.search(q, k)
        # If the chunks list and FAISS index are kept in sync
        # then we can use the returned index directly.
        # Filter out -1 indices from chunks and scores
        chunks = []
        scores = []
        for i, s in zip(Inds[0], Dists[0]):
            if i != -1:
                chunks.append(self.chunks[i])
                scores.append(s)
        return chunks, scores

    def _calculate_embeddings(self, texts: list[str]) -> list[list[float]]:
        """
        Compute embeddings for a list of texts using the configured model.

        The returned embeddings are in the same order as the input texts.

        Parameters
        ----------
        texts : list[str]
            Text strings to embed.

        Returns
        -------
        list[list[float]]
            List of embedding vectors, one per input text.

        Notes
        -----
        - All embeddings have identical dimensionality.
        - This method does not normalise embeddings; normalisation is
        performed by the caller before adding to or querying the index.
        """
        response = self.client.embeddings.create(model=self.model, input=texts)
        # Keep order of embeddings as inputs
        data = sorted(response.data, key=lambda x: x.index)
        return [item.embedding for item in data]

    def load(self, directory_path: Path | str):
        """
        Load a previously saved database from disk.

        This restores:
        - The FAISS index from 'faiss.index'
        - Metadata (including embedding model) from 'manifest.json'
        - Chunk objects from 'chunks.jsonl'

        Parameters
        ----------
        directory_path : str or Path
            Directory containing the saved database files.

        Raises
        ------
        ValueError
            If the directory does not exist or required files are missing.
            If the number of chunks does not match the number of vectors
            in the FAISS index.
        """
        directory_path = Path(directory_path)
        if not directory_path.is_dir():
            raise ValueError(f"Expected directory, got: {directory_path}")

        self.index = faiss.read_index(str(directory_path / "faiss.index"))

        self.chunks = []
        with open(directory_path / "chunks.jsonl", "r", encoding="utf-8") as f:
            for line in f:
                d = json.loads(line)
                self.chunks.append(Chunk(**d))

        manifest_path = directory_path / "manifest.json"
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)
        self.model = manifest["embedding_model"]
        # basic integrity check
        if self.index.ntotal != len(self.chunks):
            raise ValueError(
                f"Index ntotal={self.index.ntotal} != chunks={len(self.chunks)}"
            )

    def save(self, directory_path: Path | str):
        """
        Save the FAISS index, metadata, and chunks to disk.

        The following files are written to the specified directory:
        - faiss.index      : serialized FAISS index
        - manifest.json    : embedding model, dimension, and counts
        - chunks.jsonl     : one JSON-encoded Chunk per line

        Parameters
        ----------
        directory_path : str or Path
            Target directory. Created if it does not exist.

        Notes
        -----
        - This method overwrites existing files in the directory.
        - The saved database can be restored by passing the same directory
        to the constructor via the `restore` parameter.
        """
        directory_path = Path(directory_path)
        directory_path.mkdir(parents=True, exist_ok=True)

        faiss.write_index(self.index, str(directory_path / "faiss.index"))

        with open(directory_path / "chunks.jsonl", "w", encoding="utf-8") as f:
            # Write chunks on separate lines
            for ch in self.chunks:
                f.write(json.dumps(asdict(ch), ensure_ascii=False) + "\n")

        manifest = {
            "schema_version": 1,
            "embedding_model": self.model,
            "embedding_dim": int(self.index.d) if self.index is not None else None,
            "count": len(self.chunks),
            "metric": "cosine_ip_normalized",
        }
        with open(directory_path / "manifest.json", "w", encoding="utf-8") as f:
            json.dump(manifest, f, ensure_ascii=False, indent=2)
