"""MCP server for Astrocook (pilot).

Exposes a few recipes as tools of the Model Context Protocol, so that an
LLM agent (Claude Code, Claude Desktop, ...) can drive an analysis. Scope,
tools and design: ``ROADMAP.md``, section 4.1.

The server is a thin layer over the recipes: it does not change them, and
what it would need from them goes into the findings register instead.
Install it with the ``[mcp]`` extra and launch it with ``astrocook-mcp``
(``stdio`` transport).
"""
