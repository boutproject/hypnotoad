"""
LLM agent that generates tokamak meshes using Hypnotoad.
Connects to an LLM server using the openai.OpenAI API.

Usage:

    from hypnotoad import agent

    geqdsk_file = "/path/to/geqdsk.file"
    ENDPOINT_URL = "https://llm.server/v1"
    API_KEY = abdefghijklmnopqrstuvwxyz
    MODEL = "gpt-5-mini"

    hypnotoad_agent = agent.HypnotoadAgent(geqdsk_file, base_url=ENDPOINT_URL,
                                           api_key=API_KEY, model=MODEL)

In a script or command line use the `chat()` method:

    hypnotoad_agent.chat("Generate a mesh, adjusting settings until a mesh is successfully generated.")

In a Jupyter notebook the `chat_nb()` method will use `ipywidgets` and `IPython`
packages to display the output of the model and tool execution:

    hypnotoad_agent.chat_nb("Generate a mesh, adjusting settings until a mesh is successfully generated.")

Further questions and instructions can be issued:

    hypnotoad_agent.chat_nb("Please adjust settings to generate a mesh with psinorm_sol = 1.09")

Meshes generated are stored in `mesh_history`. The last mesh is available:

    hypnotoad_agent.last_mesh

To quickly plot the last mesh:

    hypnotoad_agent.plot_last_mesh()
"""

from .hypnotoad_agent import HypnotoadAgent

__all__ = ["HypnotoadAgent"]
