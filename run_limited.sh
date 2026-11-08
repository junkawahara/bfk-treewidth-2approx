#!/bin/bash
# Run a command with memory <= 32GB (virtual) and time <= 10 minutes.
ulimit -v 33554432
exec timeout 600 "$@"
