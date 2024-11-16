import os


def convert_multiline_fasta_to_oneline(
        input_fasta: str,
        output_fasta: str = None
):
    """
    Converts a multi-line FASTA file to a single-line per sequence.
    Args:
        input_fasta (str): Path to the input multi-line FASTA file.
        output_fasta (str, optional): Path to the output single-line.
            If not provided, creates a file with '_one_line.fasta' suffix.
    """
    if output_fasta is None:
        base_name, ext = os.path.splitext(input_fasta)
        output_fasta = base_name + "_one_line.fasta"
    with open(input_fasta) as polystr, open(output_fasta, "a+") as monostr:
        sequence = ""
        for line in polystr:
            line = line.strip()
            if line[0] == ">":
                if sequence != "":
                    monostr.write(f"{sequence}\n")
                sequence = ""
                monostr.write(f"{line}\n")
            else:
                sequence += line
        monostr.write(f"{sequence}")


def parse_blast_output(input_file: str, output_file: str = None):
    """
    Parses BLAST output and extracts the top match for each query sequence.
    Args:
        input_file (str): Path to the BLAST output file in txt format.
        output_file (str): Path to the output file to save the top matches,
        sorted alphabetically.
    """
    if output_file is None:
        base_name, ext = os.path.splitext(input_file)
        output_file = base_name + "_protein.txt"
    top_matches = []
    with open(input_file, "r") as infile:
        is_significant = False
        for line in infile:
            line = line.strip()
            if "Description" in line:
                is_significant = True
                continue
            if is_significant and line:
                match_description = line[:66]
                top_matches.append(match_description)
                is_significant = False
    top_matches = sorted(top_matches, key=lambda s: s.casefold())
    with open(output_file, "a+") as outfile:
        for match in top_matches:
            outfile.write(f"{match}\n")
