"""
Harmonizer registry and factory functions
"""

from .kg2 import harmonize_kg2
from .spoke import harmonize_spoke

HARMONIZERS = {
    'kg2': harmonize_kg2,
    'spoke': harmonize_spoke
}


def get_harmonizer(source_name: str):
    """Get the harmonizer function for a given source"""
    if source_name not in HARMONIZERS:
        raise ValueError(f"No harmonizer found for source: {source_name}")

    return HARMONIZERS[source_name]