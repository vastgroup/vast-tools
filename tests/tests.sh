#!/bin/bash

set -ueo pipefail

DIRBASE=$1
DIROUT=$2

# Compare with expression outcome file
cmp --silent $DIRBASE/Test-100.cRPKM $DIROUT/expr_out/Test-100.cRPKM || exit 1 

exit 0