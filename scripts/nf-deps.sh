#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <install|update>"
  exit 1
fi

action="$1"
if [[ "$action" != "install" && "$action" != "update" ]]; then
  echo "Error: Invalid action '$action'. Must be 'install' or 'update'."
  exit 1
fi

# Run the desired action for both subworkflows and modules
nf-core subworkflows \
  --git-remote git@github.com:ljwoods2/af3-nf-tools.git \
  $action af3 --dir workflows --no-preview

nf-core modules \
  --git-remote git@github.com:ljwoods2/af3-nf-tools.git \
  $action af3 --dir workflows --no-preview