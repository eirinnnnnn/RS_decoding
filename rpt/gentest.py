#!/usr/bin/env python3
"""
gen_vectors.py – build random (received, answer) pairs for your RS(63,42) project.

* Each vector:
      msg → cw = rs_encode(msg) →
      add ≤E errors + ≤T erasures (E,T obey  t0 + 2·t1 ≤ 21) →
      received line with "*" marking erasures
* Generates:
      input.txt   – 63 tokens per line (ints or "*")
      answer.txt  – <K> then the K information symbols

Change MAX_ERR, MAX_ERA on the CLI to explore the corner of the decoding radius.
"""
import argparse, random, sys
from pathlib import Path

import encoder           # rs_encode
import decoder           # rs_decode  (used only for self-check – optional)
import gf                # EXP_TABLE gives field size

# ---------- code parameters picked up from your own files --------------------
N = encoder.N            # 63
K = encoder.K            # 42
R = N - K                # 21
Q = 64                   # elements 0..63   (1 zero + 63 non-zero)

# ---------- helpers ----------------------------------------------------------
def rand_msg():
    return [random.randrange(Q) for _ in range(K)]

def add_noise(cw, max_err, max_era):
    """Return (received, erasure_positions) with t0+2t1 ≤ R."""
    while True:
        erasures = random.sample(range(N), max_era)
        # erasures = random.sample(range(N), random.randint(0, max_era))
        errs_ok  = max_err
        # errors   = [i for i in random.sample(range(N), random.randint(0, errs_ok))
        #             if i not in erasures]

        errors   = [i for i in random.sample(range(N), errs_ok)
                    if i not in erasures]
        # if 1:
        # if len(erasures) + 2*len(errors) <= R:
        if len(erasures) + 2*len(errors) == max_era + 2*max_err:
        # if len(erasures) + 2*len(errors) > R:
            break

    rx = cw.copy()
    for i in errors:
        rx[i] ^= random.randrange(1, Q)      # non-zero error
    for i in erasures:
        rx[i] = "*"                          # erasure mark for I/O

    return rx, erasures                      # (received, list[int])

# ---------- I/O helpers ------------------------------------------------------
def cw_to_line(vec):        # vec may contain ints & "*"
    return " ".join(str(x) for x in vec)

def msg_to_line(msg):
    return " ".join(map(str, msg))
    # return " ".join(map(str, [K] + msg))

# ---------- main -------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("-n", "--num", type=int, default=100,
                   help="how many test vectors")
    p.add_argument("-e", "--max_err", type=int, default=4,
                   help="max random symbol errors per vector")
    p.add_argument("-E", "--max_era", type=int, default=4,
                   help="max random erasures per vector")
    p.add_argument("-o", "--out", type=Path, default=Path("."), help="output dir")
    p.add_argument("--selfcheck", action="store_true",
                   help="verify decoder reproduces each message")
    args = p.parse_args(argv)

    first_bad = None            # <-- NEW
    in_lines, ans_lines = [], []
    for nummm in range(args.num):
        m  = rand_msg()
        cw = encoder.rs_encode(m)
        rx, eras = add_noise(cw, args.max_err, args.max_era)

        if args.selfcheck:
            # try:
            #     decoded = decoder.rs_decode(
            #         [0 if x=="*" else x for x in rx], eras)
            #     assert decoded == m
            # except Exception as exc:
            #     print("❌ self-check failed:", exc, file=sys.stderr)
            #     sys.exit(1)
            try:
                decoded = decoder.rs_decode(
                    [0 if x == "*" else x for x in rx], eras)
                assert decoded == m
            except Exception as exc:
                first_bad = (rx, eras, m, exc)
                break                    # stop after the first failure


        in_lines.append(cw_to_line(rx))
        ans_lines.append(msg_to_line(m))
    if first_bad:
        rx, eras, msg, exc = first_bad
        print("\n💥 First failing vector:\n"
            f"num{nummm}: received = {rx}\nerasures = {eras}\nmessage  = {msg}\n"
            f"decoder threw: {exc}")
        # print("Run:\n"
        #     "python3 - <<'PY'\n"
        #     "from decoder import debug_rs_decode\n"
        #     f"debug_rs_decode({[0 if x=='*' else x for x in rx]}, {eras})\n"
        #     "PY")
        return


    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "input.txt").write_text("\n".join(in_lines) + "\n")
    (args.out / "answer.txt").write_text("\n".join(ans_lines) + "\n")
    print(f"✅  Wrote {args.num} vectors to {args.out}/input.txt and answer.txt")

if __name__ == "__main__":
    main()
