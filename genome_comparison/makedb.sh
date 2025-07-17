#!/bin/bash

# May need blast to be added to path for this to work
for f in GCF*/*.faa
do
    # Strip directory and extension to create a clean base name
    base_name=$(basename "${f%.*}")

    # Create output path without quotes and with correct variable expansion
    output_path="$HOME/rotation2/genome_comparison/proteome_database/$base_name"

    # Run makeblastdb with proper output path and quoting
    makeblastdb -in "$f" -dbtype prot -title "$base_name" -out "$output_path"

    # Copy the original file into the proteome_database directory
    cp "$f" "$HOME/rotation2/genome_comparison/proteome_database/"
done