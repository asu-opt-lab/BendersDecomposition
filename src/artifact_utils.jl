# Q. It seems better to move this file to src/utils/utilsArtifact.jl.

# Artifact utilities for BendersX benchmark data
# Downloads data from GitHub Releases on first use

using Pkg.Artifacts

const ARTIFACTS_TOML = joinpath(@__DIR__, "..", "Artifacts.toml")

"""
    get_artifact_path(name::String) -> String

Download artifact if needed and return path to data directory.
"""
get_artifact_path(name::String) = ensure_artifact_installed(name, ARTIFACTS_TOML)
