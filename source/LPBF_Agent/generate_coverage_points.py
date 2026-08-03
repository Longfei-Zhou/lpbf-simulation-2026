#!/usr/bin/env python3
"""Canonical CLI wrapper for the package coverage-point generator."""

from lpbf_score.coverage_generator import main


if __name__ == "__main__":
    # main() returns statistics, not an exit code; exceptions signal failures.
    main()
