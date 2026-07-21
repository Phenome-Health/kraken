# kraken/scripts/gen_build_info_schema.py
"""Regenerate build_info.schema.json from the BuildInfo model.

Run after changing BuildInfo:  uv run python scripts/gen_build_info_schema.py
The committed schema is the contract Kestrel and biomapper2 reference.
"""

import json
from pathlib import Path

from kraken.build_info import BuildInfo

SCHEMA_PATH = Path(__file__).resolve().parents[1] / "build_info.schema.json"


def main() -> None:
    schema = BuildInfo.model_json_schema()
    SCHEMA_PATH.write_text(json.dumps(schema, indent=2) + "\n")
    print(f"Wrote {SCHEMA_PATH}")


if __name__ == "__main__":
    main()
