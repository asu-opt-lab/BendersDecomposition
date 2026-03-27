# BendersX.jl Documentation

Welcome to the documentation for **BendersX.jl**.

## Introduction
**BendersX.jl** is a modular and extensible framework for implementing Benders decomposition algorithms in Julia. The name **BendersX** reflects the design philosophy of the package:"**"X" denotes any extension beyond the basic form of Benders decomposition**, with the goal of making such extensions easy to implement, combine, and evaluate through a plug-and-play architecture.

To support this philosophy, BendersX.jl separates the algorithm into three core components—**Master**, **Oracle**, and **Environment**—each with a well-defined responsibility and interface. These components admit multiple interchangeable variants, allowing users to plug in implementations from the package library or introduce new ones without modifying the rest of the code. This plug-and-play structure supports rapid prototyping, reproducible experimentation, and fair comparison across Benders variants.

Built on top of **JuMP**, BendersX.jl allows users to formulate master and subproblem models using standard JuMP modeling syntax while delegating all algorithmic operations to the framework. The library includes a broad collection of built-in oracles (classical, unified, Pareto-optimal, split cuts, and problem-specific variants) and multiple environment (sequential, in-out stabilized, branch-and-bound). Users can easily integrate custom oracles or environments to explore new algorithmic ideas.

BendersX.jl also provides a growing suite of benchmarking examples—including facility location variants and network interdiction models—designed to support reproducible computational studies and easy comparison across algorithms.

Whether you are developing new Benders decomposition techniques, testing Benders variants, or building scalable optimization applications, BendersX.jl offers a clean, flexible, and extensible platform for working with Benders decomposition in Julia.