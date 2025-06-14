#!/usr/bin/env python3
"""
23andMe Genetic Data to TOML Converter
Reads 23andMe TSV file and converts first 1000 SNPs to TOML format
"""

import hashlib
import sys
import os
from typing import List, Dict, Tuple

def calculate_quality_metrics(markers: List[Dict]) -> Tuple[float, float, float]:
    """
    Calculate genetic data quality metrics using encoded alleles
    Returns: (call_rate, heterozygosity_rate, ti_tv_ratio)
    """
    if not markers:
        return 0.0, 0.0, 0.0
    
    total_markers = len(markers)
    missing_calls = 0
    heterozygous = 0
    transitions = 0
    transversions = 0
    
    for marker in markers:
        allele1 = marker['allele1']  # Now these are encoded numbers
        allele2 = marker['allele2']
        
        # Count missing calls (0 indicates missing data)
        if allele1 == 0 or allele2 == 0:
            missing_calls += 1
        elif allele1 in [1, 2, 3, 4] and allele2 in [1, 2, 3, 4]:
            # Count heterozygous calls
            if allele1 != allele2:
                heterozygous += 1
                
                # Count transitions vs transversions
                # Sort alleles for consistent comparison
                a1, a2 = sorted([allele1, allele2])
                
                # Transitions: A<->G (1<->3), C<->T (2<->4)
                if (a1 == 1 and a2 == 3) or (a1 == 2 and a2 == 4):
                    transitions += 1
                # Transversions: all other combinations
                elif (a1, a2) in [(1, 2), (1, 4), (2, 3), (3, 4)]:
                    transversions += 1
    
    # Calculate rates
    valid_calls = total_markers - missing_calls
    call_rate = (valid_calls / total_markers) if total_markers > 0 else 0.0
    heterozygosity_rate = (heterozygous / valid_calls) if valid_calls > 0 else 0.0
    ti_tv_ratio = (transitions / transversions) if transversions > 0 else 0.0
    
    return call_rate, heterozygosity_rate, ti_tv_ratio

def generate_challenge_hash(data: List[Dict]) -> str:
    """
    Generate a deterministic challenge hash from the encoded genetic data
    """
    # Create a string representation of the encoded data for hashing
    data_string = ""
    for marker in data:
        data_string += f"{marker['rsid']}{marker['chromosome']}{marker['position']}{marker['allele1']}{marker['allele2']}"
    
    # Generate SHA-256 hash
    hash_obj = hashlib.sha256(data_string.encode('utf-8'))
    return hash_obj.hexdigest()

def encode_allele(allele: str) -> int:
    """
    Encode allele character to number for circuit computation
    A=1, T=2, G=3, C=4, 0/-/D/I=0 (missing)
    """
    allele_encoding = {
        'A': 1, 'T': 2, 'G': 3, 'C': 4,
        '0': 0, '-': 0, 'D': 0, 'I': 0
    }
    return allele_encoding.get(allele.upper(), 0)

def encode_chromosome(chromosome: str) -> int:
    """
    Encode chromosome to number
    1-22=1-22, X=23, Y=24, MT/M=25
    """
    if chromosome.isdigit():
        chr_num = int(chromosome)
        if 1 <= chr_num <= 22:
            return chr_num
    elif chromosome.upper() == 'X':
        return 23
    elif chromosome.upper() == 'Y':
        return 24
    elif chromosome.upper() in ['MT', 'M']:
        return 25
    
    return 0  # Invalid chromosome

def validate_genetic_marker(fields: List[str]) -> bool:
    """
    Validate that a genetic marker has proper format
    """
    if len(fields) != 5:
        return False
    
    rsid, chromosome, position, allele1, allele2 = fields
    
    # Validate rsid format (should start with 'rs' or 'i')
    if not (rsid.startswith('rs') or rsid.startswith('i')):
        return False
    
    # Validate chromosome (1-22, X, Y, MT)
    if encode_chromosome(chromosome) == 0:
        return False
    
    # Validate position (should be numeric)
    try:
        pos = int(position)
        if pos <= 0:
            return False
    except ValueError:
        return False
    
    # Validate alleles (should be A, T, G, C, 0, -, D, or I)
    valid_alleles = {'A', 'T', 'G', 'C', '0', '-', 'D', 'I'}
    if allele1 not in valid_alleles or allele2 not in valid_alleles:
        return False
    
    return True

def parse_23andme_file(filename: str, max_markers: int = 1000) -> List[Dict]:
    """
    Parse 23andMe TSV file and extract genetic markers
    """
    markers = []
    
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        if not lines:
            raise ValueError("File is empty")
        
        # Skip header comments (lines starting with #)
        data_start = 0
        for i, line in enumerate(lines):
            if not line.strip().startswith('#'):
                data_start = i
                break
        
        if data_start >= len(lines):
            raise ValueError("No data found in file")
        
        # Parse header line
        header_line = lines[data_start].strip()
        headers = [h.strip() for h in header_line.split('\t')]
        
        # Validate expected headers
        expected_headers = ['rsid', 'chromosome', 'position', 'allele1', 'allele2']
        if not all(h in headers for h in expected_headers):
            raise ValueError(f"Invalid headers. Expected: {expected_headers}, Found: {headers}")
        
        print(f"✅ Valid headers found: {headers}")
        
        # Parse data lines
        processed_count = 0
        skipped_count = 0
        
        for i, line in enumerate(lines[data_start + 1:], 1):
            if processed_count >= max_markers:
                print(f"✅ Reached maximum limit of {max_markers} markers")
                break
                
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            
            fields = [f.strip() for f in line.split('\t')]
            
            # Validate marker format
            if not validate_genetic_marker(fields):
                skipped_count += 1
                if skipped_count <= 10:  # Only show first 10 warnings
                    print(f"⚠️  Skipping invalid marker on line {i}: {fields}")
                continue
            
            rsid, chromosome, position, allele1, allele2 = fields
            
            markers.append({
                'rsid': rsid,
                'chromosome': encode_chromosome(chromosome),
                'position': int(position),
                'allele1': encode_allele(allele1),
                'allele2': encode_allele(allele2),
                'raw_chromosome': chromosome,  # Keep original for reference
                'raw_allele1': allele1,        # Keep original for reference
                'raw_allele2': allele2         # Keep original for reference
            })
            
            processed_count += 1
            
            # Progress indicator
            if processed_count % 100 == 0:
                print(f"📊 Processed {processed_count} markers...")
        
        print(f"✅ Successfully processed {processed_count} markers")
        if skipped_count > 0:
            print(f"⚠️  Skipped {skipped_count} invalid markers")
        
        return markers
        
    except FileNotFoundError:
        raise FileNotFoundError(f"❌ File not found: {filename}")
    except Exception as e:
        raise Exception(f"❌ Error reading file: {str(e)}")

def write_toml_file(markers: List[Dict], output_filename: str):
    """
    Write genetic markers to TOML file with encoded alleles and chromosomes
    """
    if not markers:
        raise ValueError("No markers to write")
    
    # Calculate quality metrics
    call_rate, het_rate, ti_tv = calculate_quality_metrics(markers)
    
    # Generate challenge hash
    challenge_hash = generate_challenge_hash(markers)
    
    print(f"📊 Quality Metrics:")
    print(f"   Call Rate: {call_rate:.4f} ({call_rate*100:.2f}%)")
    print(f"   Heterozygosity Rate: {het_rate:.4f} ({het_rate*100:.2f}%)")
    print(f"   Ti/Tv Ratio: {ti_tv:.4f}")
    print(f"   Challenge Hash: {challenge_hash[:16]}...")
    
    try:
        with open(output_filename, 'w', encoding='utf-8') as file:
            # Write header comment explaining encoding
            file.write('# Genetic Data Encoding:\n')
            file.write('# Alleles: A=1, T=2, G=3, C=4, Missing=0\n')
            file.write('# Chromosomes: 1-22=1-22, X=23, Y=24, MT=25\n')
            file.write('\n')
            
            # Write header with quality metrics
            file.write(f'challenge_hash = "{challenge_hash}"\n')
            file.write(f'min_call_rate = "{call_rate:.6f}"\n')
            file.write(f'min_heterozygosity_rate = "{het_rate:.6f}"\n')
            file.write(f'ti_tv_ratio = "{ti_tv:.6f}"\n')
            file.write('\n')
            
            # Write each DNA marker with encoded values
            for i, marker in enumerate(markers):
                file.write('[[dna]]\n')
                file.write(f'allele1 = {marker["allele1"]}\n')  # Now numbers
                file.write(f'allele2 = {marker["allele2"]}\n')  # Now numbers
                file.write(f'chromosome = {marker["chromosome"]}\n')  # Now numbers
                file.write(f'position = {marker["position"]}\n')
                file.write(f'rsid = "{marker["rsid"]}"\n')
                
                # Add original values as comments for reference
                file.write(f'# Original: {marker["raw_allele1"]}/{marker["raw_allele2"]} on chr{marker["raw_chromosome"]} at {marker["position"]}\n')
                file.write('\n')
        
        print(f"✅ Successfully wrote {len(markers)} encoded markers to {output_filename}")
        print(f"🔢 Allele encoding: A=1, T=2, G=3, C=4, Missing=0")
        print(f"🧬 Chromosome encoding: 1-22=1-22, X=23, Y=24, MT=25")
        
    except Exception as e:
        raise Exception(f"❌ Error writing TOML file: {str(e)}")

def display_encoding_info():
    """
    Display the encoding scheme used for alleles and chromosomes
    """
    print("🔢 Encoding Scheme:")
    print("   Alleles: A=1, T=2, G=3, C=4, Missing/Unknown=0")
    print("   Chromosomes: 1-22=1-22, X=23, Y=24, MT=25")
    print("   This encoding optimizes genetic data for ZK circuit computation")
    print()

def validate_toml_output(filename: str):
    """
    Basic validation of the generated TOML file with encoded values
    """
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            content = file.read()
        
        # Count [[dna]] sections
        dna_sections = content.count('[[dna]]')
        
        # Check for required fields
        required_fields = ['challenge_hash', 'min_call_rate', 'min_heterozygosity_rate', 'ti_tv_ratio']
        missing_fields = [field for field in required_fields if field not in content]
        
        # Check for encoded values (should contain numbers without quotes for alleles/chromosomes)
        has_encoded_alleles = 'allele1 = 1' in content or 'allele1 = 2' in content or 'allele1 = 3' in content or 'allele1 = 4' in content
        has_encoded_chromosomes = any(f'chromosome = {i}' in content for i in range(1, 26))
        
        if missing_fields:
            print(f"⚠️  Warning: Missing fields in TOML: {missing_fields}")
        elif not has_encoded_alleles:
            print(f"⚠️  Warning: Alleles don't appear to be properly encoded as numbers")
        elif not has_encoded_chromosomes:
            print(f"⚠️  Warning: Chromosomes don't appear to be properly encoded as numbers")
        else:
            print(f"✅ TOML validation passed: {dna_sections} DNA entries with proper encoding")
        
        return dna_sections > 0 and not missing_fields and has_encoded_alleles and has_encoded_chromosomes
        
    except Exception as e:
        print(f"❌ Error validating TOML: {str(e)}")
        return False

def main():
    """
    Main function to handle command line arguments and process files
    """
    print("🧬 23andMe to TOML Converter (with Numeric Encoding)")
    print("=" * 50)
    
    # Display encoding information
    display_encoding_info()
    
    # Handle command line arguments
    if len(sys.argv) < 2:
        print("Usage: python genetic_to_toml.py <input_file.txt> [output_file.toml] [max_markers]")
        print("Example: python genetic_to_toml.py genome_data.txt output.toml 1000")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'genetic_data.toml'
    max_markers = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    
    # Validate input file exists
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file '{input_file}' does not exist")
        sys.exit(1)
    
    try:
        # Parse the genetic data file
        print(f"📖 Reading genetic data from: {input_file}")
        print(f"🎯 Processing maximum {max_markers} markers")
        markers = parse_23andme_file(input_file, max_markers)
        
        if not markers:
            print("❌ No valid genetic markers found in file")
            sys.exit(1)
        
        # Write TOML output
        print(f"💾 Writing TOML output to: {output_file}")
        write_toml_file(markers, output_file)
        
        # Validate output
        print("🔍 Validating output file...")
        if validate_toml_output(output_file):
            print("🎉 Conversion completed successfully!")
            
            # Show file sizes
            input_size = os.path.getsize(input_file)
            output_size = os.path.getsize(output_file)
            print(f"📁 Input file size: {input_size:,} bytes")
            print(f"📁 Output file size: {output_size:,} bytes")
        else:
            print("❌ Output validation failed")
            sys.exit(1)
            
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
