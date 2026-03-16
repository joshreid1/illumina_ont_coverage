import json
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python parse_summary.py <input_json> <output_tsv>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    # Load JSON data
    with open(input_path, "r") as f:
        data = json.load(f)

    # Extract required fields
    tp = data.get("TP-base", 0)
    fp = data.get("FP", 0)
    fn = data.get("FN", 0)

    # Write TSV output
    with open(output_path, "w") as out:
        out.write("Variant\tTP\tFP\tFN\n")
        out.write(f"RECORDS\t{tp}\t{fp}\t{fn}\n")

if __name__ == "__main__":
    main()
