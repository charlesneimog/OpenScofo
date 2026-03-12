import os
import csv

def main():
    base_folder = input("Enter the path of the base folder: ").strip()
    output_csv = input("Enter the CSV filename to save (e.g., output.csv): ").strip()

    results = []

    for root, dirs, files in os.walk(base_folder):
        # Check if the current folder has no subfolders
        if not dirs:
            relative_path = os.path.relpath(root, base_folder)
            print(f"Folder: {relative_path}")
            label = ""
            while label not in ["e", "o", "n"]:
                label = input("Enter label (extended=e, ordinario=o, none=n): ").strip().lower()
            results.append([relative_path, label])

    # Write results to CSV
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Folder", "Label"])
        writer.writerows(results)

    print(f"Classification saved to {output_csv}")

if __name__ == "__main__":
    main()
