import argparse
import pysam
import re
import pandas as pd
from collections import defaultdict
import sys

def process_bam_to_matrix(bam_file, output_prefix, mapq_threshold):
    print(f"Reading BAM file: {bam_file} with MAPQ threshold >= {mapq_threshold}...")
    
    # Data structure: counts[gene][(CellID, Cx, Cy)] = set(UMIs)
    expression_data = defaultdict(lambda: defaultdict(set))
    
    # Regular expression for parsing read headers
    # Example: @E...|Cx:i:14376|Cy:i:16337_5310712A3251_D07_kraken:taxid|1783272
    # Group 1: Cx | Group 2: Cy | Group 3: CellID (CID) | Group 4: UMI (MID) | Group 5: TaxID
    header_pattern = re.compile(r'Cx:i:(\d+)\|Cy:i:(\d+)_([A-Z0-9]+)_([A-Z0-9]+)_kraken:taxid\|(\d+)')
    
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Error opening BAM file: {e}")
        sys.exit(1)
        
    mapped_reads = 0
    parsed_reads = 0
    
    for read in bam.fetch():
        # Apply user-defined MAPQ threshold
        if read.is_unmapped or read.mapping_quality < mapq_threshold:
            continue
            
        mapped_reads += 1
        qname = read.query_name
        target_gene = read.reference_name
        
        match = header_pattern.search(qname)
        if match:
            cx = int(match.group(1))
            cy = int(match.group(2))
            cellid = match.group(3) # Extract the hexadecimal CID as the CellID
            umi = match.group(4)    # Extract the UMI used for deduplication
            
            # Use the combination of CellID, Cx and Cy as the unique identifier of a spatial spot
            expression_data[target_gene][(cellid, cx, cy)].add(umi)
            parsed_reads += 1

    bam.close()
    print(f"Total mapped reads passed MAPQ >= {mapq_threshold} filter: {mapped_reads}")
    print(f"Total successfully parsed with spatial info: {parsed_reads}")
    
    print("Generating spatial expression matrix...")
    records = []
    for gene, spatial_dict in expression_data.items():
        for (cellid, cx, cy), umis in spatial_dict.items():
            records.append({
                "CellID": cellid,
                "Cx": cx,
                "Cy": cy,
                "Gene": gene,
                "UMI_Count": len(umis) # Expression level after UMI deduplication
            })
            
    df = pd.DataFrame(records)
    
    # Dynamically generate the output filename, including the MAPQ threshold
    output_csv = f"{output_prefix}_MAPQ{mapq_threshold}.csv"
    
    if not df.empty:
        # Pivot table: use CellID, Cx and Cy as the row index and Gene as the columns
        matrix_df = df.pivot(index=['CellID', 'Cx', 'Cy'], columns='Gene', values='UMI_Count').fillna(0)
        # Write to file; the index columns are saved to the CSV automatically
        matrix_df.to_csv(output_csv)
        print(f"Matrix successfully saved to {output_csv}")
        print(f"Matrix shape: {matrix_df.shape[0]} spatial spots, {matrix_df.shape[1]} genes.")
    else:
        print("Warning: No data to write. Check if reads successfully mapped or if MAPQ threshold is too strict.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Stereo-seq BAM to Spatial Matrix with filtering.")
    
    # Define command-line arguments
    parser.add_argument("-i", "--input", required=True, help="Input BAM file (e.g., mapped_GUS.sorted.bam)")
    parser.add_argument("-o", "--outprefix", default="Spatial_GUS_Matrix", help="Prefix for output CSV (default: Spatial_GUS_Matrix)")
    parser.add_argument("-q", "--mapq", type=int, default=10, help="MAPQ threshold (default: 10)")
    
    args = parser.parse_args()
    
    process_bam_to_matrix(args.input, args.outprefix, args.mapq)