#!/usr/bin/env python3
"""Thin launcher; all parsing/configuration/scoring is shared with the CLI."""
import sys
from hgtfinder.cli import main

if __name__ == '__main__':
    args = sys.argv[1:]
    config = []
    if '--config' in args:
        i = args.index('--config')
        config = args[i:i+2]
        del args[i:i+2]
    main(config + ['debug-score'] + args)
