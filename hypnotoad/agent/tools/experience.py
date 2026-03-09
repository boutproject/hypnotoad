"""
Data and routines for summarizing mesh generation experience.
Intended to build a database that improves future tasks.
"""

SUMMARY_PROMPT = """
You are writing an “Experience Chunk” for a RAG knowledge base used by an LLM agent that operates the hypnotoad mesh generator.

Goal: Produce a compact, highly searchable, technically accurate summary of this run that will help future agents solve similar mesh-generation tasks. The chunk will be embedded for semantic search and also indexed for keyword/BM25 search.

Write ONLY the Experience Chunk text. Do not include JSON, code blocks, or extra commentary.

INPUTS YOU WILL RECEIVE (conceptually):
- run_status: "success" or "fail"
- equilibrium_summary: topology, number of X-points, any notable geometry/topology notes
- goal_summary: what the agent was trying to achieve (resolution goals, speed vs quality tradeoffs, etc.)
- defaults: the default options (from OptionsFactory)
- effective_options: the full options used for this run (after defaults + overrides)
- overrides: a dict of options that differ from defaults (already computed)
- diagnostics: key metrics and checks (mesh sizes, quality metrics, runtime, warnings, errors)
- artifacts: paths/URLs to config files, logs, and mesh output (if any)

STYLE AND CONTENT RULES
- Be concise but information-dense. Prefer short lines and bullets.
- Include option names exactly as they appear in the configuration (preserve nesting/paths).
- Focus on what changed from defaults, what symptoms were observed, and why the changes helped.
- Avoid speculation. If you don’t know why something helped, say “reason unclear”.
- Include keywords that improve retrieval: topology terms, common warning/error phrases, and key option names.
- If run_status is "fail", emphasize the error signature and last attempted overrides, and suggest the most plausible next changes (max 3) grounded in the observed failure.

OUTPUT FORMAT (follow exactly)
Line 1: [SUCCESS] or [FAIL] | topology=<...> | goal=<short> | eq=<short fingerprint or identifier> | version=<git sha or version if available>

Section: Situation
- 2–5 bullets describing the equilibrium/topology and the objective.

Section: Key overrides (diff from defaults)
- Group overrides by subsystem if possible (e.g., “geometry”, “spacing”, “x-point handling”, “smoothing”, “solver/integration”).
- List 5–20 overrides max. Each line:
  - <option_path>: <value>  (default: <default_value>) — <1 short clause describing intent>
- If there are more than 20 overrides, include the 20 most consequential and add one line:
  - (N more overrides omitted)

Section: Observations and lessons
- 3–8 bullets, each must follow this pattern:
  - Symptom: <what was observed>
    Change: <option(s) changed>
    Outcome: <what improved/what happened>
    Why: <brief rationale> (or “reason unclear”)
- If run_status is "fail", replace “Outcome” with “Result” and focus on:
  - error message / warning text (quote short fragments, <= 15 words)
  - where it occurred (stage: parsing, equilibrium, region detection, mesh generation, smoothing, output)
  - the most plausible next changes (max 3)

Section: Results (or Failure details)
- If success: include the most important diagnostics:
  - runtime, mesh dimensions, region count, min cell size (if available), quality checks summary
  - warnings (if any) as a short list
- If fail:
  - error_signature: <concise stable identifier; include key exception class/message fragment>
  - last_good_state: <if any> else “none”
  - what_to_try_next: 1–3 bullets (must be concrete option edits or checks)

Section: Artifacts
- config: <path or url if available>
- log: <path or url if available>
- output: <path or url if available>
- notes: <optional 1 line; e.g. “replay by running …” but keep it short>

IMPORTANT
- Do not include any sensitive or irrelevant information.
- Do not include raw stack traces; only short error fragments.
- Ensure the chunk is self-contained: a future reader should understand what worked/failed and what to try next.
"""
