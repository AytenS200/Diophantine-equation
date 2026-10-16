import argparse
import json
import sys

from . import __version__
from .core import solve
from .parser import ParseError


def _run(text, domain, limit, as_json):
    try:
        res = solve(text, domain=domain, limit=limit)
    except (ParseError, ValueError, OverflowError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    if as_json:
        print(json.dumps(res.to_dict(samples=limit), ensure_ascii=False, indent=2))
    else:
        print(res)
    return 0


def _interactive(domain, limit, as_json):
    print(f"Diophantine solver {__version__}. Enter an equation "
          "(empty line or 'q' to quit).")
    print("Commands: ':domain Z|N|P', ':limit N', ':json on|off'")
    while True:
        try:
            line = input(f"[{domain}] > ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        if line in ("", "q", "quit", "exit"):
            break
        if line.startswith(":"):
            cmd, _, arg = line[1:].partition(" ")
            arg = arg.strip()
            if cmd == "domain" and arg.upper() in ("Z", "N", "P"):
                domain = arg.upper()
            elif cmd == "limit" and arg.isdigit():
                limit = int(arg)
            elif cmd == "json" and arg in ("on", "off"):
                as_json = arg == "on"
            else:
                print("Unknown command.")
            continue
        _run(line, domain, limit, as_json)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="diophantine",
        description="Solve Diophantine equations, systems and congruences with step-by-step explanations.")
    ap.add_argument("equation", nargs="*",
                    help='equation in quotes, e.g. "3x + 5y = 7"; omit for interactive mode')
    ap.add_argument("-d", "--domain", default="Z", type=str.upper, choices=["Z", "N", "P"],
                    help="Z: all integers (default), N: non-negative, P: positive")
    ap.add_argument("-l", "--limit", type=int, default=10,
                    help="number of examples shown per infinite family (default 10)")
    ap.add_argument("--json", action="store_true", help="print the result as JSON")
    ap.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")
    a = ap.parse_args(argv)
    if a.equation:
        return _run(" ".join(a.equation), a.domain, a.limit, a.json)
    if not sys.stdin.isatty():
        code = 0
        for line in sys.stdin:
            if line.strip():
                code |= _run(line.strip(), a.domain, a.limit, a.json)
        return code
    return _interactive(a.domain, a.limit, a.json)


def _force_utf8_output():
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError, OSError):
            pass


def entry():
    _force_utf8_output()
    try:
        return main()
    except BrokenPipeError:
        import os
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 0


if __name__ == "__main__":
    sys.exit(entry())
