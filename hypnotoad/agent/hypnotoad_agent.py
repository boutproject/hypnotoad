import logging
import json
import pprint
from typing import Optional
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
4. Construct corrected settings, then call validate_settings.
   Fix all reported issues before proceeding.
5. Call run_hypnotoad with the corrected settings.

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

# Anthropic API format
TOOLS = [
    {
        "name": "validate_settings",
        "description": "Validate a settings dict before running Hypnotoad.",
        "input_schema": {
            "type": "object",
            "properties": {"settings": {"type": "object"}},
            "required": ["settings"],
        },
    },
    {
        "name": "run_hypnotoad",
        "description": "Run Hypnotoad mesh generator with a settings dict. Returns success/error and mesh metadata.",
        "input_schema": {
            "type": "object",
            "properties": {
                "settings": {
                    "type": "object",
                    "description": "Hypnotoad settings dictionary",
                }
            },
            "required": ["settings"],
        },
    },
    {
        "name": "inspect_mesh",
        "description": """Inspect mesh quality after generation. Use detail levels
progressively:
- 'summary'  (default): global pass/fail + warning list. Always call this first.
- 'standard': per-region statistics + connections. Call when summary has warnings.
- 'full':     worst-cell locations, metric tensors, interface continuity details.
              Call only to diagnose a specific problem identified at standard level.""",
        "input_schema": {
            "type": "object",
            "properties": {"detail": {"type": "string"}},
            "required": [],
        },
    },
    {
        "name": "list_meshes",
        "description": "List all meshes generated in this session with their "
        "index, pass/fail status, and settings",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "search_hypnotoad_options",
        "description": (
            "Search the Hypnotoad settings options reference. Returns the k most "
            "relevant options matching the query, each with its name, default value, "
            "type, allowed values, and description.\n\n"
            "Use this tool when:\n"
            "- You need to know the exact name of an option (e.g. 'what option "
            "controls poloidal spacing near the X-point?')\n"
            "- You need to know the default, type, or allowed values for a specific "
            "option before setting it\n"
            "- validate_settings has returned an unknown_key error and you want to "
            "find the correct option name\n"
            "- You are constructing a settings dict and want to check what options "
            "are available for a particular aspect of the mesh\n\n"
            "Do not guess option names. Always use this tool if you are unsure."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": (
                        "Natural language description of the option or behaviour you "
                        "are looking for. Can be a partial option name, a physical "
                        "concept, or a description of what you want to control. "
                        "Examples:\n"
                        "- 'X-point poloidal spacing'\n"
                        "- 'number of radial points in SOL'\n"
                        "- 'target plate resolution'\n"
                        "- 'nx_inter_sep'\n"
                        "- 'orthogonal mesh'"
                    ),
                },
                "k": {
                    "type": "integer",
                    "description": (
                        "Number of options to return. Default 4. Use a larger value "
                        "(up to 10) when exploring an unfamiliar area of the settings "
                        "space, or when the first results do not contain what you need."
                    ),
                    "default": 4,
                    "minimum": 1,
                    "maximum": 10,
                },
            },
            "required": ["query"],
        },
    },
]


def to_openai_tools(anthropic_tools: list[dict]) -> list[dict]:
    """Convert Anthropic-format tool definitions to OpenAI format."""
    return [
        {
            "type": "function",
            "function": {
                "name": t["name"],
                "description": t["description"],
                "parameters": t["input_schema"],
            },
        }
        for t in anthropic_tools
    ]


TOOLS_OPENAI = to_openai_tools(TOOLS)

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
    def __init__(self, gridfile, base_url: str = None, api_key: str = None, model=None):

        from openai import OpenAI
        from .tools.search import ChunkDatabase, extract_option_chunks

        self.gridfile = gridfile
        self.logger = logger.getChild(self.__class__.__name__)
        self.client = OpenAI(base_url=base_url, api_key=api_key)
        self.model = model

        # Check which models are available
        available_models = [model.id for model in self.client.models.list()]
        if model not in available_models:
            raise ValueError(
                f"Model {model} not available. Available models are {available_models}"
            )
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
                "signature": {"settings": dict},
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
        }

        # Index available options so that the LLM can query
        self.options_db = ChunkDatabase(extract_option_chunks(POSSIBLE_OPTIONS))

        # Store generated meshes
        self.mesh_history = []

    def chat(
        self, user_input: str, max_iterations: int = 20, task_handler=default_handler
    ) -> str:
        """
        task_handler(title, func, *args, **kwargs) : function
            Wrapper that should print the title and then run func(*args, **kwargs)
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
                    "tool_calls": [
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
                    else None,
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
                        def run_task():
                            print(
                                f"Calling {tc.function.name}\nInputs: {pprint.pformat(args)}",
                                flush=True,
                            )

                            try:
                                tool = self.tool_registry[tc.function.name]
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
        """Wrapper that handles model output in a Jupyter notebook"""

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

    def run_hypnotoad(self, settings: dict = {}) -> dict:
        """Run Hypnotoad with given settings dict"""
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
            self.mesh_history.append({"settings": settings, "mesh": mesh})
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
        Inspect a previously generated mesh by index.
        mesh_index: index from run_hypnotoad result. -1 = most recent (default).
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
        """Summarise all mesh attempts in this session."""
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
        Describe the magnetic equilibrium and geometry. Call this at the start
        of a session or when you need to understand the physics before choosing
        mesh settings. Returns metrics that inform resolution and spacing choices.
        """
        try:
            return tools.describe_equilibrium(self.gridfile)
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def search_hypnotoad_options(self, query: str, k: int = 4) -> list[dict]:
        """ """
        return self.options_db.retrieve(query, k)

    @property
    def last_mesh(self) -> Optional[BoutMesh]:
        """The last successfully generated mesh. Can be None."""
        if len(self.mesh_history) == 0:
            return None
        return self.mesh_history[-1]["mesh"]

    @property
    def last_settings(self) -> Optional[dict]:
        """Return the settings used to create the most recent successful mesh.
        Can be None."""
        if len(self.mesh_history) == 0:
            return None
        return self.mesh_history[-1]["settings"]

    def plot_last_mesh(self, ax=None):
        """Plots the most recent successfully generated mesh"""
        mesh = self.last_mesh
        if mesh is None:
            return
        ax = mesh.plotPotential(axis=ax)
        return mesh.plotGridCellEdges(ax=ax)
