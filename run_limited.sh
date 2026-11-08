#!/bin/bash
# SPDX-License-Identifier: MIT
# Run a command with memory <= 32GB (virtual) and time <= 10 minutes:
# SIGTERM after 590 s and, if the command ignores it, SIGKILL 10 s later.
# GNU timeout is preferred: the uutils "timeout" (coreutils 0.2.2, the default "timeout"
# on Ubuntu 25.10) did not stop the command with -k or -s KILL in our tests.
ulimit -v 33554432
TIMEOUT=$(command -v gnutimeout || command -v timeout)
exec "$TIMEOUT" -k 10 590 "$@"
