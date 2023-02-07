import sys, subprocess
from os import path, listdir


def create_skmer_fa(outdir):
    """ Merge all the txt files that contains skmers into a single fasta file called all_skmers.fa
    """
    with open(path.join(outdir, "all_skmers.fa"), 'w') as fa:
        # Read each minimizer file one by one
        for filename in listdir(outdir):
            if not filename.endswith(".txt"):
                continue

            filepath = path.join(outdir, filename)
            minimizer = filename[:-4]
            skmer_idx=0

            with open(filepath) as sk_fp:
                # Read each superkmer one by one
                for line in sk_fp:
                    line = line.strip()

                    if len(line) == 0:
                        continue

                    print(f">{minimizer}-{skmer_idx}\n{line}", file=fa)
                    skmer_idx += 1


def count(fasta_file, k):
    with open("/dev/null", "w") as nullfp:
        name = fasta_file[:-3]
        print("--- Counting kmers ---")
        subprocess.run(f"kmc -k{k} -fm -ci0 {fasta_file} {name} /tmp".split())
        print("--- Dumping kmers ---")
        subprocess.run(f"kmc_dump -ci0 {name} {name}.counts".split())
        print("--- Sorting kmers ---")
        subprocess.run(f"sort {name}.counts -o {name}.sorted".split())
        print("--- Cleaning ---")
        subprocess.run(f"rm {name}.kmc_pre {name}.kmc_suf {name}.counts".split())


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("python3 verif.py genome.fa skmer_outdir/ k")

    fa = sys.argv[1]
    outdir = sys.argv[2]
    k = int(sys.argv[3])
    sk_fa = path.join(outdir, "all_skmers.fa")

    create_skmer_fa(outdir)
    count(fa, k)
    count(sk_fa, k)

    print()
    print("--- diff origin/skmers ---")
    subprocess.run(f"diff {fa[:-3]}.sorted {sk_fa[:-3]}.sorted".split())
    subprocess.run(f"rm {fa[:-3]}.sorted {sk_fa[:-3]}.sorted".split())
