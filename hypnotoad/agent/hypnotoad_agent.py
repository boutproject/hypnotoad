import logging
import json
import pprint
from typing import Optional
from pathlib import Path
from ..cases import tokamak
from ..core.mesh import BoutMesh
from . import tools

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """
You are an expert assistant for generating 2D tokamak plasma simulation meshes 
using Hypnotoad. You help users configure mesh settings, run the mesh generator, 
diagnose problems, and iteratively improve mesh quality. The user is a plasma
physicist and you can use technical terminology freely.

## Tools Available

- get_equilibrium_info: Describe the magnetic equilibrium and geometry.
  Call this at the start of a session or when you need to understand the
  physics before choosing mesh settings. Returns metrics that inform
  psinorm, resolution and spacing choices.

- validate_settings(settings): Check a settings dict for missing required keys 
  and type errors before running.

- run_hypnotoad(settings): Run the mesh generator. Returns success/error and 
  mesh diagnostics. On error, the message will include the exception message.

- inspect_mesh(detail): Inspect the most recently generated mesh file.

- search_hypnotoad_options(query, k): Search the Hypnotoad settings options
  reference, returning the k most relevant results. Do not guess settings keys.

- search_experience(query, k): Search past experience, returning the k most
  relevant results.

## Workflow

Follow these steps in order. Do not skip steps or change their sequence.

### Step 1 — Establish initial settings

IF the user provides settings or options:
  - Call validate_settings with those settings
  - Fix every issue reported before proceeding
  - Report what was fixed to the user
ELSE:
  - Use settings = {
        'nx_core': 18,
        'nx_sol': 18,
        'orthogonal': True,
        'psinorm_pf': 0.96,
        'psinorm_sol': 1.05,
        'target_all_poloidal_spacing_length': 0.05,
        'xpoint_poloidal_spacing_length': 0.15,
    }
  - call search_hypnotoad_options to identify
    appropriate changes based on the geometry and user request.
  - Call validate_settings on the constructed settings and fix any issues

Do not call run_hypnotoad until validate_settings returns valid=true.

### Step 2 — Generate a mesh

Call run_hypnotoad with the most recent validated settings.
Use the same input to run_hypnotoad as was used to validate_settings. 

IF run_hypnotoad succeeds: go to Step 3.
IF run_hypnotoad fails:    go to Step 2a.

#### Step 2a — Diagnose and fix a failed run

Do not guess at settings changes. Follow this sequence:

1. Read the error message carefully. If it names a specific option or
   parameter, call search_hypnotoad_options with that name first.
   The error "Cannot create connected double-null grid" indicates that
   nx_inter_sep should be set to 1 or larger. 
2. Call get_equilibrium_info if you have not already done so. Pay attention to:
   - The number and psinorm values of the X-points.
     Only X-points within the normalised psi range [psinorm_pf, psinorm_sol] will
     be included in the mesh. A maximum of two X-points can be included.
     If more than two X-points fall within the domain, narrow psinorm_sol to
     exclude the excess.
   - The suggested psinorm_pf and psinorm_sol values.
3. Call search_hypnotoad_options to find options relevant to the error.
4. Before changing more than 2 options, call search_experience for similar cases.
5. Construct corrected settings, then call validate_settings.
   Fix all reported issues before proceeding.
6. Call run_hypnotoad with the corrected settings.

Repeat Step 2a up to 3 times. If the mesh still fails after 3 attempts,
report the full error history to the user and ask for guidance.
Do not continue attempting without user input.

### Step 3 — Assess mesh quality

Call inspect_mesh(detail="summary") immediately after every successful
run_hypnotoad, even if the run looked clean.

IF valid=false (negative Jacobian or self-crossing cells):
  - This mesh must not be used. Treat as a failure.
  - Call inspect_mesh(detail="standard") to identify which regions failed.
  - Go to Step 2a to fix the settings.

IF valid=true but warnings are present:
  - Call inspect_mesh(detail="standard") to investigate.
  - Report each warning to the user with a plain-language explanation of
    its physical significance.
  - Suggest specific settings changes that would address each warning.
    Support each suggestion with a search_hypnotoad_options call —
    do not suggest option names from memory.
  - Ask the user whether to proceed with improvements or accept the mesh.

IF valid=true and no warnings:
  - Report the mesh quality summary to the user.
  - Ask the user whether the mesh is acceptable or further refinement
    is needed.

    
### General rules

- Never call run_hypnotoad with settings that have not passed validate_settings.
- Never suggest or use an option name that has not been confirmed by
  search_hypnotoad_options or returned by validate_settings.
- Never make more than one settings change at a time without explaining
  the reason for each change.
- If uncertain about any step, ask the user before proceeding.
"""

TOOLS_OPENAI = [
    {
        "type": "function",
        "function": {
            "name": "get_equilibrium_info",
            "description": (
                "Describe the magnetic equilibrium and geometry from the input grid file. "
                "Call this at the start of a session (before choosing mesh settings) or when you need "
                "to understand topology (single-null / double-null), X-point locations, and size/shape "
                "metrics that inform resolution and spacing choices."
            ),
            "parameters": {"type": "object", "properties": {}, "required": []},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "validate_settings",
            "description": (
                "Validate a Hypnotoad settings dict against the OptionsFactory schema (types, allowed values, constraints). "
                "Call this BEFORE run_hypnotoad when you have changed settings or are unsure about option names/values. "
                "If validation reports unknown keys, use search_hypnotoad_options to find the correct option names."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "settings": {
                        "type": "object",
                        "description": "Hypnotoad settings dictionary to validate. Keys must be exact option paths; values must match types/constraints.",
                    },
                },
                "required": ["settings"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_hypnotoad",
            "description": (
                "Run the Hypnotoad mesh generator with a settings dict. "
                "Returns success/error and a mesh_index for later inspection. "
                "Best practice: validate_settings -> run_hypnotoad -> inspect_mesh(detail='summary') "
                "and only increase detail if needed."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "settings": {
                        "type": "object",
                        "description": "Hypnotoad settings dictionary. Use validate_settings first; do not guess option names.",
                    },
                },
                "required": ["settings"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "inspect_mesh",
            "description": (
                "Inspect mesh quality after generation.\n"
                "Use detail levels progressively:\n"
                "- 'summary'  (default): global pass/fail + warning list. Always call this first.\n"
                "- 'standard': per-region statistics + connections. Call when summary has warnings.\n"
                "- 'full':     worst-cell locations, metric tensors, interface continuity details.\n"
                "             Call only to diagnose a specific problem identified at standard level."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "mesh_index": {
                        "type": "integer",
                        "description": (
                            "Index of mesh to inspect, as returned by run_hypnotoad. "
                            "Use -1 to inspect the most recent mesh."
                        ),
                        "default": -1,
                    },
                    "detail": {
                        "type": "string",
                        "description": "Inspection detail level.",
                        "enum": ["summary", "standard", "full"],
                        "default": "summary",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "list_meshes",
            "description": (
                "List all meshes generated in this session with their mesh_index, pass/fail status, warning count, and settings."
            ),
            "parameters": {"type": "object", "properties": {}, "required": []},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_hypnotoad_options",
            "description": (
                "Search the Hypnotoad settings options reference (BM25 keyword search over option docs). "
                "Returns the k most relevant options matching the query, each with its name, default value, type, allowed values, and description.\n\n"
                "Use this tool when:\n"
                "- You need the exact name/path of an option\n"
                "- You need default/type/allowed values before setting it\n"
                "- validate_settings reports unknown keys\n"
                "- You are exploring how to control a specific aspect of the mesh\n\n"
                "Do not guess option names. Always use this tool if you are unsure."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": (
                            "Natural language description of the option/behaviour. Can be a partial option name or concept. Examples:\n"
                            "- 'X-point poloidal spacing'\n"
                            "- 'number of radial points in SOL'\n"
                            "- 'target plate resolution'\n"
                            "- 'nx_intersep'\n"
                            "- 'orthogonal mesh'"
                        ),
                    },
                    "k": {
                        "type": "integer",
                        "description": "Number of options to return. Default 4. Use up to 10 when exploring.",
                        "default": 4,
                        "minimum": 1,
                        "maximum": 10,
                    },
                },
                "required": ["query"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_experience",
            "description": (
                "Search the saved experience database of prior Hypnotoad runs (successful and failed). "
                "Use this BEFORE making large settings changes, especially when you see a warning/error or when working with a similar topology.\n\n"
                "Typical uses:\n"
                "- 'connected double-null second X-point distortion'\n"
                "- 'nx_intersep too low warnings'\n"
                "- 'mesh smoothing interface continuity'\n"
                "- paste a short error/warning message fragment to find prior fixes\n\n"
                "Returns the top-k most relevant experience reports with key overrides and lessons."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Natural-language query, option names, topology keywords, or error/warning fragments.",
                    },
                    "k": {
                        "type": "integer",
                        "description": "Number of experience records to return. Default 4. Use up to 10 when exploring.",
                        "default": 4,
                        "minimum": 1,
                        "maximum": 10,
                    },
                },
                "required": ["query"],
            },
        },
    },
]

POSSIBLE_OPTIONS = (
    tokamak.TokamakEquilibrium.user_options_factory.defaults
    | tokamak.TokamakEquilibrium.nonorthogonal_options_factory.defaults
    | BoutMesh.user_options_factory.defaults
)


def normalise_arguments(signature: dict, args):
    """Return a dict matching the given signature or raise a ValueError.

    A common issue with LLM tools is inconsistency in the argument
    format. Sometimes the arguments are passed in order, other times
    as a dict with argument names. This function attempts to normalise
    the arguments to match the signature.
    """
    # Check first argument
    if signature == {}:
        return {}  # No arguments
    first_key, first_type = next(iter(signature.items()))

    if isinstance(args, dict):
        if isinstance(args, first_type):
            # This dict could either be intended to be the first argument,
            # or to contain function arguments.
            if any([key not in signature for key in args]):
                # Assign arguments to the first key
                return {first_key: args}
        # Filter keys to those in the signature
        return {key: value for key, value in args.items() if key in signature}

    elif isinstance(args, first_type):
        return {first_key: args}
    raise ValueError(f"Arguments {args} do not match signature {signature}")


def default_handler(title, func, *args, **kwargs):
    """Wrapper for tool calls"""
    print(title)
    return func(*args, **kwargs)


class HypnotoadAgent:
    """
    LLM-driven controller for the Hypnotoad mesh generator.

    The agent exposes tools that the language model may call (validate_settings,
    run_hypnotoad, inspect_mesh, list_meshes, search_hypnotoad_options,
    search_experience, get_equilibrium_info). It maintains an in-memory
    session (messages + mesh_history) and optional persistent experience
    storage (ChunkFaissDatabase).

    Parameters
    ----------
    gridfile : str or Path
        Path to the equilibrium/grid file (GEQDSK or compatible format) used
        to construct meshes.

    base_url : str, optional
        Optional base URL for the OpenAI-compatible API.

    api_key : str, optional
        API key for the OpenAI-compatible client.

    model : str, optional
        Chat/completion model id used for agent reasoning and tools.

    embedding_model : str, optional
        Embedding model id used to create experience embeddings.

    experience_db : str or Path, optional
        Directory to restore/save the persistent experience FAISS store.

    Notes
    -----
    - The agent is primarily a thin orchestration layer; heavy lifting is
      delegated to tools in `tools.*` and to the OpenAI client for LLM calls.
    - Mesh generation results are stored in `self.mesh_history` as entries
      containing at least {'settings', 'mesh', 'diagnostics'}.
    """

    def __init__(
        self,
        gridfile,
        base_url: str = None,
        api_key: str = None,
        model: str = None,
        embedding_model: str = None,
        experience_db: Path | str = None,
    ):
        """
        Initialize the HypnotoadAgent.

        Sets up the OpenAI client, model/tool bindings, BM25 options index,
        optional FAISS-based experience database, and an empty mesh history.

        Parameters
        ----------
        gridfile : str | Path
            Path to the equilibrium/grid file for mesh generation.
        base_url : str, optional
            Base URL for the OpenAI-compatible API.
        api_key : str, optional
            API key for the OpenAI-compatible API.
        model : str, optional
            Model id for chat completions.
        embedding_model : str, optional
            Embedding model id for the experience database.
        experience_db : str | Path, optional
            Directory to restore the experience database from.
        """
        from openai import OpenAI
        from .tools.search import ChunkDatabase, extract_option_chunks

        self.gridfile = gridfile
        self.logger = logger.getChild(self.__class__.__name__)
        self.client = OpenAI(base_url=base_url, api_key=api_key)
        self.model = model
        self.tools = TOOLS_OPENAI

        # Maintain chat history. This is sent to the LLM at each call
        self.messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
        ]

        self.tool_registry = {
            "get_equilibrium_info": {
                "function": self.get_equilibrium_info,
                "signature": {},
            },
            "validate_settings": {
                "function": lambda **args: tools.validate_settings(
                    POSSIBLE_OPTIONS, **args
                ),
                "signature": {"settings": dict},
            },
            "run_hypnotoad": {
                "function": self.run_hypnotoad,
                "signature": {"settings": dict, "notes": str},
            },
            "inspect_mesh": {
                "function": self._inspect_mesh,
                "signature": {"mesh_index": int, "detail": str},
            },
            "list_meshes": {"function": self.list_meshes, "signature": {}},
            "search_hypnotoad_options": {
                "function": self.search_hypnotoad_options,
                "signature": {"query": str, "k": int},
            },
            "search_experience": {
                "function": self.search_experience,
                "signature": {"query": str, "k": int},
            },
        }

        # Index available options so that the LLM can query
        self.options_db = ChunkDatabase(extract_option_chunks(POSSIBLE_OPTIONS))

        # Database of past experience
        self.experience_db = None
        self.experience_db_path = experience_db
        if embedding_model or experience_db:
            self._init_experience_db(
                embedding_model=embedding_model, experience_db=experience_db
            )

        # Store generated meshes
        self.mesh_history = []

    def _init_experience_db(
        self,
        embedding_model: str = None,
        experience_db: Path | str = None,
    ):
        """Initialise the experience database.

        embedding_model : str, optional
            Embedding model id for the experience database.
        experience_db : str | Path, optional
            Directory to restore the experience database from.
        """
        from .tools.search import ChunkFaissDatabase

        if experience_db:
            experience_db = Path(experience_db)
            if not experience_db.is_dir():
                self.logger.warning(
                    f"Experience DB '{experience_db}' does not exist. Will be created on save."
                )
                experience_db = None  # Don't try to restore

        self.experience_db = ChunkFaissDatabase(
            self.client, model=embedding_model, restore=experience_db
        )

    def chat(
        self, user_input: str, max_iterations: int = 20, task_handler=default_handler
    ) -> str:
        """
        Drive an interactive agent loop with the LLM, handling tool calls.

        This appends the user's input to the internal message history, sends the
        conversation to the LLM, and executes any tool calls returned by the LLM.
        Tool executions are wrapped and dispatched through `task_handler` so the
        caller can capture, display, or redirect output.

        Parameters
        ----------
        user_input : str
            Natural-language instruction or question for the agent.
        max_iterations : int, default=20
            Maximum number of LLM iterations (tool-call cycles) to perform.
        task_handler : callable
            Signature: task_handler(title: str, func: Callable[[], Any]) -> Any.
            Used to run tool calls; allows UI integration (e.g., capturing output).

        Returns
        -------
        str
            The assistant's final textual reply (may be empty string).

        Notes
        -----
        - Tool call arguments are expected to be JSON strings and will be parsed.
        - Tool results are appended to the conversation as 'tool' messages so the
          LLM can continue reasoning with tool outputs.
        - This method mutates `self.messages`. Consider cloning if you want an
          ephemeral reasoning call that doesn't alter session history.
        - The function protects against malformed tool arguments but tool errors
          are returned as structured error objects to the LLM.
        """
        self.logger.debug(f"User input: {user_input}")
        self.messages.append({"role": "user", "content": user_input})

        for it in range(max_iterations):
            response = self.client.chat.completions.create(
                model=self.model,
                tools=TOOLS_OPENAI,
                messages=self.messages,
            )

            self.logger.debug(response)

            choice = response.choices[0]
            msg = choice.message

            # Add assistant message -- must keep tool_calls intact
            self.messages.append(
                {
                    "role": "assistant",
                    "content": msg.content,  # may be None
                    "tool_calls": (
                        [
                            {
                                "id": tc.id,
                                "type": "function",
                                "function": {
                                    "name": tc.function.name,
                                    "arguments": tc.function.arguments,  # keep as string
                                },
                            }
                            for tc in msg.tool_calls
                        ]
                        if msg.tool_calls
                        else None
                    ),
                }
            )

            if choice.finish_reason == "tool_calls":
                # Add tool results -- one per tool call, matched by tool_call_id
                for tc in msg.tool_calls:
                    try:
                        args = json.loads(tc.function.arguments)
                    except json.JSONDecodeError as e:
                        args = {}
                        result = {
                            "status": "error",
                            "error": "invalid_json",
                            "message": f"Tool arguments were not valid JSON: {e}",
                            "raw": tc.function.arguments,
                        }
                    else:
                        # Wrap the tool call in a function to pass to task_handler
                        # This enables output to be captured and redirected in
                        # the user interface.

                        # Use default arguments to avoid potential late-binding
                        # closure bug if task_handler defers tasks.
                        tc_name = tc.function.name

                        def run_task(tc_name=tc_name, args=args):
                            print(
                                f"Calling {tc_name}\nInputs: {pprint.pformat(args)}",
                                flush=True,
                            )

                            try:
                                tool = self.tool_registry[tc_name]
                                # Normalise the arguments to match signature
                                result = tool["function"](
                                    **normalise_arguments(tool["signature"], args)
                                )
                            except Exception as e:
                                result = {
                                    "status": "error",
                                    "message": str(e),
                                }
                            print(f"Result: {pprint.pformat(result)}", flush=True)
                            return result

                        result = task_handler(f"Calling {tc.function.name}", run_task)

                    self.messages.append(
                        {
                            "role": "tool",
                            "tool_call_id": tc.id,  # must match the id in the assistant message
                            "content": json.dumps(result),
                        }
                    )
            else:
                return msg.content or ""
        return "Exceeded maximum iterations. See log for details."

    def chat_nb(self, user_input: str, max_iterations: int = 20):
        """
        Notebook-friendly wrapper around `chat` that captures tool output in
        collapsible UI widgets (ipywidgets).

        Parameters
        ----------
        user_input : str
            User instruction to pass to the agent.
        max_iterations : int, default=20
            Maximum number of LLM iterations.

        Returns
        -------
        None
            Prints the final assistant text and presents interactive UI elements
            for tool execution logs.

        Notes
        -----
        - This method requires Jupyter/IPython (ipywidgets). It is a convenience
          wrapper and does not change agent semantics.
        """

        import ipywidgets as widgets
        from IPython.display import display

        def run_step(title, func, *args, **kwargs):
            """Runs a function and captures all output into a collapsible Accordion."""
            # Create an output widget to capture prints/logs
            out = widgets.Output()

            # Create the Accordion UI
            accordion = widgets.Accordion(children=[out], selected_index=None)
            accordion.set_title(0, f"▶ {title}")
            display(accordion)

            # Capture the output
            with out:
                # Redirect logging to this specific output widget
                handler = logging.StreamHandler()
                logger.addHandler(handler)
                try:
                    return func(*args, **kwargs)
                finally:
                    logger.removeHandler(handler)

        result = self.chat(
            user_input, task_handler=run_step, max_iterations=max_iterations
        )
        print(result)

    def run_hypnotoad(
        self, settings: Optional[dict] = None, notes: Optional[str] = None
    ) -> dict:
        """
        Run the Hypnotoad mesh generator using the provided settings.

        This method:
        - Loads the equilibrium from self.gridfile with provided settings,
        - Constructs a BoutMesh, runs the standard processing (calculateRZ,
            geometry, etc.),
        - Computes diagnostics via tools.inspect_mesh(detail='summary'),
        - Appends a dictionary to `self.mesh_history` with keys:
            {'settings', 'mesh', 'diagnostics'}.

        Parameters
        ----------
        settings : dict, optional
            Hypnotoad settings dictionary. If None, defaults are used.

        Returns
        -------
        dict
            Structured result with at minimum:
            - status: 'success' or 'error'
            - mesh_index: integer index into mesh_history (when success)
            - n_meshes: total number of saved meshes
            - diagnostics: diagnostics dict (when success)
            - message: error message (when failure)
            - hint: optional next-step hint

        Notes
        -----
        - Call validate_settings before run_hypnotoad when possible.
        - Exceptions during reading or mesh generation are caught and returned
          as structured errors (status='error').
        """
        settings = settings or {}
        if notes:
            print(notes)
        try:
            # Read the grid file
            with open(self.gridfile, "rt") as fh:
                eq = tokamak.read_geqdsk(fh, settings=settings)

            if isinstance(eq, tuple):
                raise eq[1]  # Second element is an exception
            mesh = BoutMesh(eq, settings)
            mesh.calculateRZ()
            mesh.geometry()
            idx = len(self.mesh_history)
            diagnostics = tools.inspect_mesh(mesh, detail="summary")
            self.mesh_history.append(
                {"settings": settings, "mesh": mesh, "diagnostics": diagnostics}
            )
            return {
                "status": "success",
                "mesh_index": idx,  # <-- LLM uses this for inspect_mesh
                "n_meshes": len(self.mesh_history),
                "diagnostics": diagnostics,
                "hint": f"Use inspect_mesh(mesh_index={idx}, detail='standard') to get more detail, "
                f"or inspect_mesh(mesh_index=N) for any previous mesh.",
            }
        except Exception as e:
            return {
                "status": "error",
                "mesh_index": None,
                "message": str(e),
            }

    def _inspect_mesh(self, mesh_index: int = -1, detail: str = "summary") -> dict:
        """
        Inspect a stored mesh by index and return diagnostics.

        This is a thin wrapper around tools.inspect_mesh that selects the mesh
        from `self.mesh_history`.

        Parameters
        ----------
        mesh_index : int, default=-1
            Index of the mesh to inspect. -1 selects the most recent mesh.
        detail : str, default='summary'
            Level of inspection: 'summary', 'standard', or 'full'.

        Returns
        -------
        dict
            The same structure returned by tools.inspect_mesh, or an error object
            with keys:
            - status: 'error'
            - message: error text
            - hint: optional usage hint

        Raises
        ------
        None
            All exceptions are captured and returned as structured error dicts.
        """
        if len(self.mesh_history) == 0:
            return {
                "status": "error",
                "message": "No mesh generated. Use run_hypnotoad to generate a mesh.",
            }
        try:
            mesh = self.mesh_history[mesh_index]["mesh"]
        except Exception as e:
            return {
                "status": "error",
                "message": str(e),
                "hint": f"mesh_index={mesh_index} must be in range 0 <= mesh_index < {len(self.mesh_history)}.",
            }
        return tools.inspect_mesh(mesh, detail=detail)

    def list_meshes(self) -> dict:
        """
        Return a summary of all meshes generated in this session.

        The returned object contains:
        - n_meshes: int
        - meshes: list of dicts, each containing:
            - mesh_index: int
            - valid: bool (diagnostics.get('valid', False))
            - n_warnings: int (diagnostics.get('n_warnings', 0))
            - settings: dict (the settings used to produce the mesh)

        Returns
        -------
        dict
            Session-level mesh summary.

        Notes
        -----
        - This is a lightweight listing intended for quick inspection by the LLM.
        - For in-depth diagnostics call inspect_mesh on a specific mesh_index.
        """
        return {
            "n_meshes": len(self.mesh_history),
            "meshes": [
                {
                    "mesh_index": idx,
                    "valid": a["diagnostics"].get("valid", False),
                    "n_warnings": a["diagnostics"].get("n_warnings", 0),
                    "settings": a["settings"],
                }
                for idx, a in enumerate(self.mesh_history)
            ],
        }

    def get_equilibrium_info(self) -> dict:
        """
        Describe the magnetic equilibrium and geometry associated with self.gridfile.

        This function returns physics-informed metrics that guide mesh choices,
        such as topology (single-null/double-null), X-point locations, device
        extents, and shape proxies. It is intended to be called at session start
        or before choosing mesh settings.

        Returns
        -------
        dict
            Either a description dict (topology, key coordinates, scalar metrics),
            or an error object: {'status': 'error', 'message': str}.

        Notes
        -----
        - Implementation calls tools.describe_equilibrium(self.gridfile).
        - The returned structure should be concise (a few scalars + short textual
          indicators) so that it fits well into the model context.
        """
        try:
            return tools.describe_equilibrium(self.gridfile)
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def search_hypnotoad_options(self, query: str, k: int = 4) -> list[dict]:
        """
        Search the options reference for matching Hypnotoad settings.

        This wraps the BM25-based `self.options_db` lookup and returns
        JSON-serializable option descriptors suitable for LLM consumption.

        Parameters
        ----------
        query : str
            Natural language or partial option name to search for.
        k : int, default=4
            Number of results to return.

        Returns
        -------
        list[dict]
            List of option summaries. Each dict should contain at least:
            - name/path (exact configuration key)
            - default value
            - type
            - allowed values or constraints (if known)
            - short description or example

        Notes
        -----
        - The LLM should call this before guessing option names or setting unknown keys.
        - This method returns structured dicts (not Chunk objects) to keep tool
          results easy to parse by the LLM.
        """
        return self.options_db.retrieve(query, k)

    def search_experience(self, query: str, k: int = 4):
        if self.experience_db is None:
            return []
        chunks, scores = self.experience_db.retrieve(query, k)
        return [
            {"text": c.text, "score": s, "source": c.source, "section": c.section}
            for c, s in zip(chunks, scores)
        ]

    @property
    def last_mesh(self) -> Optional[BoutMesh]:
        """
        The most recent successfully generated BoutMesh, or None.

        Returns
        -------
        BoutMesh or None
            The mesh object for programmatic inspection/plotting.
        """
        if len(self.mesh_history) == 0:
            return None
        return self.mesh_history[-1]["mesh"]

    @property
    def last_settings(self) -> Optional[dict]:
        """
        The settings dict used to generate the most recent successful mesh,
        or None if no successful mesh exists.

        Returns
        -------
        dict or None
            The resolved settings dict (defaults applied) for the last mesh.
        """
        if len(self.mesh_history) == 0:
            return None
        return self.mesh_history[-1]["settings"]

    def plot_last_mesh(self, ax=None):
        """
        Plot the most recent successfully generated mesh.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Optional axis to draw into. If None, the mesh's default plotting
            behavior will create or return an axis.

        Returns
        -------
        matplotlib.axes.Axes or None
            The axis containing the plotted mesh, or None if no mesh exists.

        Notes
        -----
        - This convenience method delegates to the BoutMesh plotting helpers:
          mesh.plotPotential() and mesh.plotGridCellEdges().
        """
        mesh = self.last_mesh
        if mesh is None:
            return
        ax = mesh.plotPotential(axis=ax)
        return mesh.plotGridCellEdges(ax=ax)

    def add_experience_report(self, embedding_model: Optional[str] = None):
        """
        Summarize the most recent run and add an 'experience' Chunk to the
        experience database.

        Behavior:
          - Constructs a compact summary (using an ephemeral LLM call) that
            includes: topology, goal, overrides (diff from defaults),
            3-6 lessons (symptom→change→outcome), and diagnostics summary.
          - Creates a Chunk(section='experience', chunk_type='experience')
            with the summary text.
          - Computes embedding(s) and adds them to the experience DB via
            self.experience_db.add_chunks([chunk]).

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If no experience DB is configured (self.experience_db is None)
            and no embedding_model is provided.

        Notes
        -----
        - Should use a one-shot LLM call (not append to self.messages) to avoid
          corrupting the ongoing conversational history.
        - The helper should compute `overrides` as the diff between the last
          settings and OptionsFactory defaults for compactness and reproducibility.
        """
        from .tools import experience
        from .tools.search import Chunk

        if self.experience_db is None:
            if embedding_model is None:
                raise ValueError(
                    "No experience DB configured and no embedding_model provided."
                )
            from .tools.search import ChunkFaissDatabase

            self.experience_db = ChunkFaissDatabase(self.client, model=embedding_model)

        # Generate a summary including key lessons learned
        summary = self.chat(experience.SUMMARY_PROMPT)

        ch = Chunk(
            text=summary,
            section="experience",
            source="experience",
            chunk_type="experience",
        )
        # Add chunk to the database
        self.experience_db.add_chunks([ch])

    def save_experience(self, path: Path | str = None):
        """
        Persist the experience database to disk.

        Parameters
        ----------
        path : str or Path, optional
            Destination directory. If None, uses the path provided at
            initialization (self.experience_db_path). If that is also None,
            a ValueError is raised.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If no destination path is provided and no experience DB path was
            configured during initialization.
        """
        if self.experience_db is None:
            return
        if path is None:
            # Use the path given to init (may be None)
            path = self.experience_db_path
        if path is None:
            raise ValueError("No path given to save_experience()")
        self.experience_db.save(path)
