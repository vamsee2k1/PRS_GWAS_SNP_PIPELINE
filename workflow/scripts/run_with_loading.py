#!/usr/bin/env python3
import argparse
import queue
import subprocess
import sys
import threading
import time


def parse_args():
    p = argparse.ArgumentParser(
        description="Run a command and show a terminal loading spinner between output bursts."
    )
    p.add_argument(
        "--interval",
        type=float,
        default=0.2,
        help="Spinner refresh interval in seconds (default: 0.2)",
    )
    p.add_argument(
        "command",
        nargs=argparse.REMAINDER,
        help="Command to run (prefix with -- before the command).",
    )
    args = p.parse_args()
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        p.error("No command provided. Example: run_with_loading.py -- snakemake -n")
    return args


def pump_output(stream, q):
    try:
        for line in iter(stream.readline, ""):
            q.put(line)
    finally:
        stream.close()


def clear_spinner(width=100):
    sys.stdout.write("\r" + (" " * width) + "\r")
    sys.stdout.flush()


def main():
    args = parse_args()
    cmd = args.command

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    q = queue.Queue()
    t = threading.Thread(target=pump_output, args=(proc.stdout, q), daemon=True)
    t.start()

    spinner = "|/-\\"
    spin_i = 0
    start = time.time()
    last_was_spinner = False

    while True:
        if proc.poll() is not None and q.empty():
            break

        try:
            line = q.get(timeout=args.interval)
            if last_was_spinner:
                clear_spinner()
                last_was_spinner = False
            sys.stdout.write(line)
            sys.stdout.flush()
        except queue.Empty:
            elapsed = int(time.time() - start)
            mins = elapsed // 60
            secs = elapsed % 60
            msg = f"[{spinner[spin_i % len(spinner)]}] working... {mins:02d}:{secs:02d}"
            sys.stdout.write("\r" + msg)
            sys.stdout.flush()
            spin_i += 1
            last_was_spinner = True

    if last_was_spinner:
        clear_spinner()
    t.join(timeout=0.2)
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
