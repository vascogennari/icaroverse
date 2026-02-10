"""
Use a job template file
"""
import os
import argparse


def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)


def main():
    parser = argparse.ArgumentParser()
    input_group = parser.add_argument_group(
        "Input"
    )
    input_group.add_argument(
        '-t', '--template',
        required=True,
        type=argparse.FileType(mode="r"),
        help="job template file",
    )
    mandatory_group = parser.add_argument_group(
        "Mandatory substitution"
    )
    mandatory_group.add_argument(
        '-e', '--executable',
        required=True,
        help="executable to use"
    )
    mandatory_group.add_argument(
        '-c', '--config',
        required=True,
        help="config file to use"
    )
    output_group = parser.add_argument_group(
        "Output"
    )
    output_group.add_argument(
        '-o', '--output',
        required=True,
        type=argparse.FileType(mode="w+"),
        help="job file",
    )
    output_group.add_argument(
        '-s', '--submit',
        choices=['slurm', 'condor', 'none'],
        default="none",
        help="submit job (use none to no submit)"
    )

    args, extra_args = parser.parse_known_args()
    # 1. Read template
    template = args.template.read()
    # 2. Create substitution dictionary
    # Mandatory arguments
    d = {
        "executable": args.executable,
        "config":     args.config,
    }
    # Extra-arguments
    for k, v in pairwise(extra_args):
        k = k.lstrip(parser.prefix_chars)
        d[k] = v
    # 3. Substitute
    output = template.format(**d)
    # 4. Write output
    args.output.writelines(output)
    # Submit job
    if args.submit == 'slurm':
        os.system(f"sbatch {args.output}")
    if args.submit == 'condor':
        os.system(f"condor_submit {args.output}")


if __name__ == '__main__':
    main()
