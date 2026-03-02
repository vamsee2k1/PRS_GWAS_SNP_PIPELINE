#!/usr/bin/env python3
import argparse
import queue
import re
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
        "--clean-snakemake",
        action="store_true",
        help="Show concise Snakemake progress in terminal (full details remain in Snakemake logs).",
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


class SnakemakeCleanFormatter:
    TIMESTAMP_RE = re.compile(r"^\[[A-Za-z]{3}\s+[A-Za-z]{3}\s+\d+\s+\d{2}:\d{2}:\d{2}\s+\d{4}\]\s*$")
    RULE_START_RE = re.compile(r"^(?:localrule|rule)\s+([A-Za-z0-9_]+):\s*$")
    RULE_DONE_RE = re.compile(r"^Finished jobid:\s+\d+\s+\(Rule:\s+([A-Za-z0-9_]+)\)\s*$")
    STEP_RE = re.compile(r"^\d+\s+of\s+\d+\s+steps\s+\(\d+%\)\s+done\s*$")
    EXECUTE_RE = re.compile(r"^Execute\s+\d+\s+jobs\.\.\.\s*$")
    ERROR_RE = re.compile(
        r"^(Error in rule|RuleException:|WorkflowError:|CalledProcessError|Cannot convert input table|Exiting because a job execution failed)"
    )
    TOTAL_RE = re.compile(r"^\s*total\s+(\d+)\s*$")

    def __init__(self):
        self.pending_timestamp = ""

    def _consume_timestamp(self):
        ts = self.pending_timestamp
        self.pending_timestamp = ""
        return ts

    def transform(self, line):
        txt = line.rstrip("\n")
        if not txt:
            return None

        if self.TIMESTAMP_RE.match(txt):
            self.pending_timestamp = txt
            return None

        m = self.RULE_START_RE.match(txt)
        if m:
            ts = self._consume_timestamp()
            rule = m.group(1)
            prefix = f"{ts} " if ts else ""
            return f"{prefix}START {rule}\n"

        m = self.RULE_DONE_RE.match(txt)
        if m:
            ts = self._consume_timestamp()
            rule = m.group(1)
            prefix = f"{ts} " if ts else ""
            return f"{prefix}DONE {rule}\n"

        if self.STEP_RE.match(txt):
            self.pending_timestamp = ""
            return f"{txt}\n"

        if self.EXECUTE_RE.match(txt):
            self.pending_timestamp = ""
            return f"{txt}\n"

        m = self.TOTAL_RE.match(txt)
        if m:
            self.pending_timestamp = ""
            return f"Total jobs: {m.group(1)}\n"

        if txt.startswith("Config file "):
            self.pending_timestamp = ""
            return f"{txt}\n"
        if txt.startswith("Assuming unrestricted shared filesystem usage."):
            self.pending_timestamp = ""
            return f"{txt}\n"
        if txt.startswith("host: "):
            self.pending_timestamp = ""
            return f"{txt}\n"
        if txt.startswith("Building DAG of jobs..."):
            self.pending_timestamp = ""
            return f"{txt}\n"
        if txt.startswith("Using shell: "):
            self.pending_timestamp = ""
            return f"{txt}\n"
        if txt.startswith("Provided cores: "):
            self.pending_timestamp = ""
            return f"{txt}\n"
        if txt.startswith("Complete log(s): "):
            self.pending_timestamp = ""
            return f"{txt}\n"

        if self.ERROR_RE.search(txt):
            self.pending_timestamp = ""
            return f"{txt}\n"
        if txt.startswith("[ERROR]"):
            self.pending_timestamp = ""
            return f"{txt}\n"

        return None


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
    formatter = SnakemakeCleanFormatter() if args.clean_snakemake else None

    spinner = "|/-\\"
    spin_i = 0
    start = time.time()
    last_was_spinner = False

    while True:
        if proc.poll() is not None and q.empty():
            break

        try:
            line = q.get(timeout=args.interval)
            if formatter is not None:
                line = formatter.transform(line)
                if line is None:
                    continue
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
