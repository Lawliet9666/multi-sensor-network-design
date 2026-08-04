# Contributing

Thank you for helping improve this paper-reproduction repository.

1. Open an issue describing the result, script, or documentation to change.
2. Keep changes focused on the seven paper experiments and the CPU Julia path.
3. Run `julia --project=. -e 'using Pkg; Pkg.test()'`.
4. Run `julia --project=. scripts/reproduce_paper.jl --profile smoke`.
5. In the pull request, report both commands and the generated artifact path.

Please do not add CUDA, Python notebooks, unrelated kernels, or large generated
results without first discussing the scope. Generated runs belong under
`results/`.
