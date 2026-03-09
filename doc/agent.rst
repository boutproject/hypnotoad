HypnotoadAgent: LLM-Assisted Mesh Generation
============================================

The ``HypnotoadAgent`` provides an LLM-driven interface to the Hypnotoad
mesh generator. It combines:

- Structured tool calls (validate, run, inspect)
- A searchable options reference
- A growing database of past experience
- Interactive or notebook-based workflows

The agent is designed to help users explore the settings space efficiently,
diagnose mesh quality issues, and improve results over time.

Overview
--------

The agent wraps Hypnotoad’s normal workflow:

1. Inspect equilibrium geometry
2. Choose mesh settings
3. Validate settings
4. Run mesh generation
5. Inspect mesh quality
6. Iterate until acceptable

The LLM can call tools to perform these steps programmatically, rather than
guessing option names or values.

Basic Usage
-----------

Create an agent by providing a GEQDSK (or equivalent) equilibrium file:

.. code-block:: python

   from hypnotoad.agent import HypnotoadAgent

   agent = HypnotoadAgent(
       gridfile="example.geqdsk",
       model="gpt-4o-mini",
       embedding_model="text-embedding-3-large",
       experience_db="experience_store"
   )

Then interact with it:

.. code-block:: python

   response = agent.chat(
       "Generate a mesh for this equilibrium with high resolution near the X-point."
   )
   print(response)

In a Jupyter notebook, use:

.. code-block:: python

   agent.chat_nb("Create a double-null mesh with refined SOL resolution")

This displays tool calls and outputs in collapsible sections.

Main Tools
----------

The agent exposes the following tools to the LLM:

get_equilibrium_info
^^^^^^^^^^^^^^^^^^^^

Describes magnetic topology, X-point locations, and geometric scale.
Call this at the start of a session to inform resolution choices.

validate_settings
^^^^^^^^^^^^^^^^^

Checks a settings dictionary against Hypnotoad’s options schema:

- Unknown option names
- Type mismatches
- Violations of allowed ranges
- Constraint failures

Always validate before running Hypnotoad.

run_hypnotoad
^^^^^^^^^^^^^

Runs the mesh generator with a given settings dictionary.

Returns:

- ``status`` (success or error)
- ``mesh_index`` (for later inspection)
- ``diagnostics`` (summary metrics and warnings)

inspect_mesh
^^^^^^^^^^^^

Examines mesh quality at different detail levels:

- ``summary`` — global pass/fail and warnings
- ``standard`` — per-region statistics and connectivity
- ``full`` — cell-level diagnostics

Typical workflow:

.. code-block:: python

   validate_settings(...)
   run_hypnotoad(...)
   inspect_mesh(detail="summary")

search_hypnotoad_options
^^^^^^^^^^^^^^^^^^^^^^^^^

Searches the settings reference using BM25 keyword search.

Use this when:

- You do not know the exact option name
- Validation reports an unknown key
- You want to understand default values or allowed ranges

Do not guess option names.

search_experience
^^^^^^^^^^^^^^^^^

Searches previously stored mesh-generation experience.

Use this when:

- You encounter a warning or error
- You are working with a similar topology
- You want to see what worked before

This enables the agent to improve over time.

Experience Database
-------------------

The agent can store summaries of successful (and failed) runs in a
persistent FAISS-based vector database.

After generating a mesh, call:

.. code-block:: python

   agent.add_experience_report()
   agent.save_experience()

Each experience entry includes:

- Topology and goal
- Overrides (differences from defaults)
- Key lessons (symptom → change → outcome)
- Diagnostics summary

On future runs, the agent can retrieve similar experiences and reuse
successful parameter combinations.

Session History
---------------

All meshes generated in a session are stored in memory:

.. code-block:: python

   agent.list_meshes()

You can inspect or plot the most recent mesh:

.. code-block:: python

   agent.inspect_mesh(mesh_index=-1, detail="standard")
   agent.plot_last_mesh()

Best Practices
--------------

1. Always call ``get_equilibrium_info`` at the start.
2. Use ``search_hypnotoad_options`` rather than guessing option names.
3. Validate settings before running.
4. Inspect at ``summary`` level before requesting ``full`` detail.
5. Change only a few options per iteration.
6. Save experience after successful runs.

Limitations
-----------

- The agent does not automatically delete meshes from history.
- The experience database is append-only.
- Large conversations may exceed model context limits.
- Mesh quality metrics depend on Hypnotoad’s diagnostics.

Future Extensions
-----------------

Possible enhancements include:

- Hybrid search (BM25 + vector) over options and experience
- Automatic option-diff computation for experience entries
- Objective-driven optimisation loops
- Integration with CI pipelines

