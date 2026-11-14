#!/usr/bin/env bash
set -euo pipefail
bash analysis/01_mafft_and_consensus.sh
python analysis/02_numbering_audit_4B7R.py
