from typing import List
import pathlib
import time
import subprocess
import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", type=int, required=True)
    parser.add_argument("-m", type=int, required=True)
    parser.add_argument("-o", "--out_dir", type=str, default="out")
    parser.add_argument(
        "-r",
        "--erlb",
        type=str,
        default="0.0",
        help="Encoding rate lower bound",
    )
    parser.add_argument("--dlb", type=int, default=0, help="Distance lower bound")
    parser.add_argument(
        "-n",
        "--nis",
        type=int,
        default=10000,
        help="Number of information sets used for distance upper bound estimation",
    )
    parser.add_argument(
        "-c",
        "--chunk_size",
        type=int,
        default=5000,
        help="Chunk size for each searching subtask",
    )
    parser.add_argument(
        "-t",
        "--target_k",
        type=int,
        default=0,
        help="Target number of logical qubits for searching",
    )
    parser.add_argument(
        "-p",
        type=int,
        default=1,
        help="Number of logical processors used for hpc-gap",
    )
    parser.add_argument(
        "-f",
        "--format_csv",
        action="store_true",
        help="Format the code coefficients in the csv file",
    )
    args = parser.parse_args()
    erlb = eval(args.erlb)

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)
    out_csv_filepath = (
        out_dir
        / f"l{args.l}_m{args.m}_erlb{erlb:.2f}_dlb{args.dlb}_nis{args.nis}_chunk{args.chunk_size}_targetk{args.target_k}.csv"
    )
    gap_call_scripts = f"""l := {args.l};
m := {args.m};
resultsFilepath := "{out_csv_filepath}";
encodingRateLowerBound := {erlb};
distanceLowerBound := {args.dlb};
numInformationSets := {args.nis};
chunkSize := {args.chunk_size};
targetK := {args.target_k};
SearchBBCodes(l, m, resultsFilepath, encodingRateLowerBound, distanceLowerBound, numInformationSets, chunkSize, targetK);
QUIT;
"""
    # print("Running GAP with the following script:")
    # print(gap_call_scripts)
    start_time = time.monotonic()
    subprocess.run(
        [
            "gap",
            "-P",
            f"{args.p}",
            "--quitonbreak",
            "-b",
            "-q",
            "search.g",
            "-c",
            f"{gap_call_scripts}",
        ]
    )
    elapsed_time = time.monotonic() - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")

    if args.format_csv:
        with open(out_csv_filepath, "r", newline="") as f:
            csv.reader(f)
            header = next(f).strip().split(",")
            rows = [row for row in csv.reader(f)]
            formatted_rows = [
                [
                    _format_coefficients(eval(row[0])),
                    _format_coefficients(eval(row[1])),
                    row[2],
                    eval(row[3]),
                ]
                for row in rows
            ]
        with open(out_csv_filepath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(formatted_rows)


def _format_coefficients(A: List[List[int]]) -> str:
    return "+".join([f'{['x', 'y'][i]}^{p}' for i, p in A])


if __name__ == "__main__":
    main()
