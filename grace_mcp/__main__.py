"""Run the GRACE MCP server: python -m grace_mcp --data-dir /path/to/datasets"""
from grace_mcp.server import main

if __name__ == "__main__":   # guard: training workers re-import this module
    main()
