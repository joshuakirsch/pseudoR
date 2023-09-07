#!/usr/bin/env python3

import argparse
import tempfile
from concurrent.futures import ThreadPoolExecutor
from shutil import which
import subprocess
import threading
import sys
import re
import os

def run_prodigal(opts, workDir, currentId, chunkFile):
    cmd = ["prodigal", "-q", "-i", chunkFile.name ]

    if opts.proteins:
        cmd.append("-a")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".faa")

    if opts.closed:
        cmd.append("-c")

    if opts.nucl:
        cmd.append("-d")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".fna")

    if opts.format:
        cmd.append("-f")
        cmd.append(opts.format)

    if opts.gencode:
        cmd.append("-g")
        cmd.append(str(opts.gencode))

    if opts.mask:
        cmd.append("-m")

    if opts.nosd:
        cmd.append("-n")

    if opts.procedure:
        cmd.append("-p")
        cmd.append(opts.procedure)

    if opts.scorefile:
        cmd.append("-s")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".score")

    cmd.append("-o")
    cmd.append(workDir.name + "/chunk" + str(currentId) + ".out")

    #print(str(cmd))
    subprocess.run(cmd, shell=False, check=True)

def append_fasta_file(file, startNum, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                if line[0] == '>':
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4) + "\n"
                trgt.write(line)
    return startNum

def append_gff_file(file, startNum, seqnumStart, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    seqnumpattern = re.compile(r"(.*)seqnum=(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                line = line.rstrip("\n")
                if line[0] == '#' and "seqnum=" in line:
                    match = re.match(seqnumpattern, line)
                    seqnumStart = seqnumStart + 1
                    line = match.group(1) + "seqnum=" + str(seqnumStart) + match.group(3)

                elif line[0] != '#' and "ID=" in line:
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4)

                # print gff header only for first chunk
                if seqnumStart != 0 and "gff-version" in line:
                    continue

                trgt.write(line + "\n")
    return (startNum, seqnumStart)


def append_gbk_file(file, startNum, seqnumStart, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    seqnumpattern = re.compile(r"(.*)seqnum=(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                line = line.rstrip("\n")

                # renumber DEFINITON line
                if line[0] == 'D' and "seqnum=" in line:
                    match = re.match(seqnumpattern, line)
                    seqnumStart = seqnumStart + 1
                    line = match.group(1) + "seqnum=" + str(seqnumStart) + match.group(3)

                elif line[0] == ' ' and "ID=" in line:
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4)

                trgt.write(line + "\n")

    return (startNum, seqnumStart)


def append_raw_file(file, targetFile, seqnumStart):
    pattern = re.compile(r"(.*)seqnum=(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                line = line.rstrip("\n")
                if "seqnum=" in line:
                    match = re.match(pattern, line)
                    seqnumStart = seqnumStart + 1
                    line = match.group(1) + "seqnum=" + str(seqnumStart) + match.group(3)
                trgt.write(line + "\n")
    return seqnumStart



def print_gff_file(file, startNum):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(file, "r") as input:
        for line in input:
            if line[0] != '#' and "ID=" in line:
                match = re.match(pattern, line)
                if match and match.group(3) == "1":
                    startNum = startNum + 1
                line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4)
            print(line)
    return startNum


def print_gbk_file(file, startNum, seqnumStart):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    seqnumpattern = re.compile(r"(.*)seqnum=(\d+)(.*)")
    with open(file, "r") as input:
        for line in input:
            line = line.rstrip("\n")

            # renumber DEFINITON line
            if line[0] == 'D' and "seqnum=" in line:
                match = re.match(seqnumpattern, line)
                seqnumStart = seqnumStart + 1
                line = match.group(1) + "seqnum=" + str(seqnumStart) + match.group(3)

            elif line[0] == ' ' and "ID=" in line:
                match = re.match(pattern, line)
                if match and match.group(3) == "1":
                    startNum = startNum + 1
                line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4)

            print(line)

    return (startNum, seqnumStart)



def print_raw_file(file, seqnumStart):
    pattern = re.compile(r"(.*)seqnum=(\d+)(.*)")
    with open(file, "r") as input:
        for line in input:
            line = line.rstrip("\n")
            if "seqnum=" in line:
                match = re.match(pattern, line)
                seqnumStart = seqnumStart + 1
                line = match.group(1) + "seqnum=" + str(seqnumStart) + match.group(3)
            print(line)
    return seqnumStart


def main():
    argp=argparse.ArgumentParser(description='Parallel Prodigal gene prediction')
    argp.add_argument('-a', "--proteins", type=str, help="Write protein translations to the selected file.")
    argp.add_argument('-c', "--closed", action="store_true", help="Closed ends.  Do not allow genes to run off edges.")
    argp.add_argument('-d', "--nucl", type=str, help="Write nucleotide sequences of genes to the selected file.")
    argp.add_argument('-f', "--format", type=str, help="Select output format (gbk, gff, or sco).  Default is gbk.")
    argp.add_argument('-g', "--gencode", type=int, help="Specify a translation table to use (default 11).")
    argp.add_argument('-i', "--input", type=str, help="Specify FASTA/Genbank input file (default reads from stdin).")
    argp.add_argument('-m', "--mask", action="store_true", help="Treat runs of N as masked sequence; don't build genes across them.")
    argp.add_argument('-n', "--nosd", action="store_true", help="Bypass Shine-Dalgarno trainer and force a full motif scan.")
    argp.add_argument('-o', "--output", type=str, help="Specify output file (default writes to stdout).")
    argp.add_argument('-p', "--procedure", type=str, help="Select procedure (single or meta).  Default is single.")
    argp.add_argument('-s', "--scorefile", type=str, help="Write all potential genes (with scores) to the selected file.")
    argp.add_argument('-T', "--tasks", type=int, help="number of prodigal processes to start in parallel (default: 20)")
    argp.add_argument('-C', "--chunksize", type=int, help="number of input sequences to process within a chunk (default: 2000)")
    opts = argp.parse_args()

    # if invoked without arguments, display usage
    if len(sys.argv) == 1:
        argp.print_help()
        exit(0)

    tasks = 20
    if opts.tasks is not None:
        if opts.tasks < 1:
            raise ValueError
        tasks = opts.tasks

    if which("prodigal") is None:
        raise ValueError("prodigal not found!")

    if opts.chunksize and opts.chunksize < 1:
        raise ValueError

    seqsPerChunk = 2000
    if opts.chunksize is not None:
        if opts.chunksize < 1:
            raise ValueError
        seqsPerChunk = opts.chunksize

    seqCnt = 0
    currentChunk = 1

    workDir = tempfile.TemporaryDirectory()
    executor = ThreadPoolExecutor(max_workers=tasks)
    currentFile = open(workDir.name + "/chunk" + str(currentChunk), 'w')

    queryFile = None
    if not opts.input:
        queryFile = "/dev/fd/0"
        if sys.stdin.isatty():
            print("Cannot read sequences from STDIN.")
            exit(1)
    else:
        queryFile = opts.input


    with open(queryFile, 'r') as fasta:
        for line in fasta:

            if line[0] == '>' and seqCnt == seqsPerChunk:
                currentFile.close()
                executor.submit(run_prodigal, opts, workDir, currentChunk, currentFile)
                currentFile = None
                seqCnt = 0
                currentChunk = currentChunk + 1

            if currentFile is None:
                currentFile = open(workDir.name + "/chunk" + str(currentChunk), 'w')

            currentFile.write(line)

            if line[0] == '>':
                seqCnt = seqCnt + 1

    if seqCnt > 0:
        currentFile.close()
        executor.submit(run_prodigal, opts, workDir, currentChunk, currentFile)

    # await completion of tasks
    executor.shutdown(wait=True)

    # collect output
    #
    proteinFile = opts.proteins
    nuclFile = opts.nucl
    outFile = opts.output
    scoreFile = opts.scorefile


    # remove output files
    if proteinFile and os.path.isfile(proteinFile):
        os.remove(proteinFile)
    if nuclFile and os.path.isfile(nuclFile):
        os.remove(nuclFile)
    if outFile and os.path.isfile(outFile):
        os.remove(outFile)
    if scoreFile and os.path.isfile(scoreFile):
        os.remove(scoreFile)


    protIdStart = 0
    nuclIdStart = 0
    gffIdStart = 0
    gbkIdStart = 0
    seqnumStart = 0
    for cur in range(1, currentChunk + 1):
        if proteinFile:
            protIdStart = append_fasta_file(workDir.name + "/chunk" + str(cur) + ".faa", protIdStart, proteinFile)
        if nuclFile:
            nuclIdStart = append_fasta_file(workDir.name + "/chunk" + str(cur) + ".fna", nuclIdStart, nuclFile)
        if scoreFile:
            seqnumStart = append_raw_file(workDir.name + "/chunk" + str(cur) + ".score", scoreFile, seqnumStart)

        if outFile:
            if opts.format == "gff":
                (gffIdStart, seqnumStart) = append_gff_file(workDir.name + "/chunk" + str(cur) + ".out", gffIdStart, seqnumStart, outFile)
            elif opts.format == "sco":
                seqnumStart = append_raw_file(workDir.name + "/chunk" + str(cur) + ".out", outFile, seqnumStart)
            else:
                (gbkIdStart, seqnumStart) = append_gbk_file(workDir.name + "/chunk" + str(cur) + ".out", gbkIdStart, seqnumStart, outFile)
        else:
            if opts.format == "gff":
                gffIdStart = print_gff_file(workDir.name + "/chunk" + str(cur) + ".out", gffIdStart)
            elif opts.format == "sco":
                seqnumStart = print_raw_file(workDir.name + "/chunk" + str(cur) + ".out", seqnumStart)
            else:
                (gbkIdStart, seqnumStart) = print_gbk_file(workDir.name + "/chunk" + str(cur) + ".out", gbkIdStart, seqnumStart)


if __name__== "__main__":
    main()
