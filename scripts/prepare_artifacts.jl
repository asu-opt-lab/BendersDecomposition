#!/usr/bin/env julia
# Script to prepare dataset artifacts for BendersX
# Usage: julia --project=. scripts/prepare_artifacts.jl
#
# This script:
# 1. Creates tar.gz archives of each dataset directory
# 2. Computes git-tree-sha1 and sha256 hashes
# 3. Outputs updated Artifacts.toml entries

using Pkg
using Pkg.Artifacts
using SHA
using Tar
using CodecZlib

const PROBLEMS_DIR = joinpath(@__DIR__, "..", "src", "problems")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "artifacts_output")
const ARTIFACTS_TOML = joinpath(@__DIR__, "..", "Artifacts.toml")

# Dataset definitions: (artifact_name, problem_type, subdirectory)
const DATASETS = [
    ("cflp_random_data", "cflp", "random_data"),
    ("cflp_locssall", "cflp", "locssall"),
    ("cflp_output", "cflp", "output"),
    ("uflp_locssall", "uflp", "locssall"),
    ("uflp_allkoerkelghosh", "uflp", "AllKoerkelGhosh"),
    ("scflp", "scflp", "SCFLP"),
    ("snip", "snip", "SNIP"),
    ("mcndp", "mcndp", "NDR"),
]

function create_artifact_archive(name::String, source_dir::String, output_dir::String)
    if !isdir(source_dir)
        @warn "Source directory not found, skipping: $source_dir"
        return nothing
    end
    
    mkpath(output_dir)
    tarball_path = joinpath(output_dir, "$(name).tar.gz")
    
    # Create tarball
    @info "Creating archive for $name..."
    open(tarball_path, "w") do io
        gzip = GzipCompressorStream(io)
        Tar.create(source_dir, gzip)
        close(gzip)
    end
    
    # Compute sha256
    sha256_hash = open(tarball_path, "r") do io
        bytes2hex(sha256(io))
    end
    
    # Compute git-tree-sha1
    # Note: Pkg.GitTools is not exported by default, we access it via Pkg
    tree_hash = Pkg.GitTools.tree_hash(source_dir)
    git_tree_sha1 = bytes2hex(tree_hash)
    
    return (
        name = name,
        tarball = tarball_path,
        sha256 = sha256_hash,
        git_tree_sha1 = git_tree_sha1,
    )
end

function main()
    println("=" ^ 60)
    println("BendersX Artifact Preparation Script")
    println("=" ^ 60)
    println()
    
    results = []
    
    for (name, problem_type, subdir) in DATASETS
        source_dir = joinpath(PROBLEMS_DIR, problem_type, "data", subdir)
        result = create_artifact_archive(name, source_dir, OUTPUT_DIR)
        if result !== nothing
            push!(results, result)
        end
    end
    
    println()
    println("=" ^ 60)
    println("SUMMARY")
    println("=" ^ 60)
    println()
    println("Generated $(length(results)) artifact archives in: $OUTPUT_DIR")
    println()
    
    # Print Artifacts.toml entries
    println("Update the following entries in Artifacts.toml:")
    println("-" ^ 60)
    
    for r in results
        println("""
[$( r.name)]
git-tree-sha1 = "$(r.git_tree_sha1)"
lazy = true

    [[$(r.name).download]]
    url = "https://github.com/asu-opt-lab/BendersX.jl/releases/download/data-v1/$(r.name).tar.gz"
    sha256 = "$(r.sha256)"
""")
    end
    
    println("-" ^ 60)
    println()
    println("Next steps:")
    println("1. Upload the .tar.gz files from $OUTPUT_DIR to GitHub Releases")
    println("2. Update the download URLs in Artifacts.toml if using a different host")
    println("3. Copy the git-tree-sha1 and sha256 values to Artifacts.toml")
    println()
    
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
