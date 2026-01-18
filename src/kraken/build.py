#!/usr/bin/env python3
"""
KRAKEN Build System
Main entry point for building the unified knowledge graph
"""

from kraken.orchestrator import KrakenBuildOrchestrator


def main():
    KrakenBuildOrchestrator().run()


if __name__ == "__main__":
    main()
