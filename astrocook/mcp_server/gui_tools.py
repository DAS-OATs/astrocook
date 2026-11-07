"""The GUI launched by the server (phase 6-ter, PR 5-bis; decision 7 of the
phase, option b): ``gui_open`` saves a session to an archive and launches the
Astrocook GUI on it; the store saves the archive again after every change of
the session, and the GUI reloads it, so that the user sees the analysis as it
goes. The session is read only in the GUI, which can make an editable copy.

What the GUI and the server share (the archive, its marker, whether the GUI
can run, the command) is in :mod:`astrocook.follow`; the tool is a thin
wrapper.
"""
from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError

from astrocook import follow
from astrocook.mcp_server.outputs import GuiResult
from astrocook.mcp_server.session_tools import SESSION_ID
from astrocook.mcp_server.sessions import SessionStore


def register_gui_tools(server: MCPServer, store: SessionStore) -> None:
    """Add the tools of the GUI to ``server``."""

    @server.tool()
    def gui_open(session_id: SESSION_ID) -> GuiResult:
        """Show a session to the user in the Astrocook GUI, on the machine
        where the server runs, and keep it up to date: after every change
        of the session (a recipe, an undo or redo, a note) the GUI shows the
        new state, keeping the user's view. The user sees it read only (they
        can make an editable copy in the GUI); nothing they do there reaches
        you. Sessions opened later by a recipe are not followed: call
        gui_open on them if the user should see them. Call it again if the
        user closed the GUI. Where there is no display (e.g. the Docker
        image) it says that the GUI is not available."""
        reason = follow.gui_unavailable_reason()
        if reason is not None:
            raise ToolError(f"The GUI is not available here: {reason} Save the "
                            "session with session_save and open the archive in "
                            "Astrocook on a machine with a display.")
        with store.lock:
            store.get(session_id)
            try:
                path = store.follow(session_id)
            except (OSError, ValueError) as e:
                raise ToolError(f"The archive for the GUI could not be written: {e}") from e
        try:
            command = follow.launch_gui(path)
        except OSError as e:
            raise ToolError(f"The GUI could not be launched: {e}. The archive "
                            f"is {path}: it can be opened in Astrocook.") from e
        return GuiResult(session_id=session_id, path=str(path), command=command)
