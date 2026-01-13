# Shared artifact utilities for BendersLibrary
# This module provides a unified way to access data through Julia Artifacts

using Pkg.Artifacts

# Path to Artifacts.toml
const ARTIFACTS_TOML = joinpath(@__DIR__, "..", "Artifacts.toml")

"""
    get_artifact_path(artifact_name::String, fallback_subdir::String, problem_dir::String) -> String

Get the path to artifact data. Uses local data if available, otherwise downloads from artifact.

# Arguments
- `artifact_name`: Name of the artifact in Artifacts.toml (e.g., "snip")
- `fallback_subdir`: Subdirectory name under data/ (e.g., "SNIP")  
- `problem_dir`: The @__DIR__ of the calling module

# Returns
Path to the data directory (either local or downloaded artifact)
"""
function get_artifact_path(artifact_name::String, fallback_subdir::String, problem_dir::String)
    # First try local path
    local_path = joinpath(problem_dir, "data", fallback_subdir)
    if isdir(local_path)
        return local_path
    end
    
    # Try to use artifact
    if isfile(ARTIFACTS_TOML)
        try
            # Use correct API: ensure_artifact_installed(name, artifacts_toml_path)
            artifact_path = ensure_artifact_installed(artifact_name, ARTIFACTS_TOML)
            return artifact_path
        catch e
            @warn "Failed to load artifact $artifact_name: $e"
        end
    end
    
    error("Data not found for $artifact_name. Either restore local data or ensure Artifacts.toml is configured correctly.")
end
