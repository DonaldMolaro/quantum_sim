# Debugging Guide

Quick checks and common mistakes for students.

## Common mistakes

- `Bit order confusion`:
  In this simulator, `q0` is LSB. Bitstrings are little-endian in gate indexing.
- `Target out of range`:
  Ensure target indices are in `[0, 2^n - 1]`.
- `Wrong angle units`:
  Rotations default to radians; use `DEG` suffix/token when needed.
- `Too many Grover rounds`:
  Success can drop after overshooting optimal iteration count.

## Use `CHECK` in CLI

`CHECK` helps validate state assumptions quickly:

```text
CHECK NORMALIZED
CHECK BASIS <index>
CHECK TARGETS <t1> [t2 ...]
CHECK BELL <PHI+|PHI-|PSI+|PSI->
```

Examples:

```text
CHECK NORMALIZED
CHECK TARGETS 3 5
CHECK BELL PHI+
```

If `CHECK BELL` fails, verify qubit order and phase-gate placement first.
