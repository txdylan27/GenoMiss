"""Evidence modules, one per line. Each exposes run(hit, cfg, log) -> dict result
and degrades gracefully (returns status='skipped' with a reason) when its optional
inputs are missing."""
