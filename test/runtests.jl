using Test
using MPSToolkit

@testset "module loads" begin
  @test isdefined(MPSToolkit, :scarfinder_step!)
end

include("test_core.jl")
include("test_docstrings.jl")
include("test_finite_tebd.jl")
include("test_finite_tdvp.jl")
include("test_operator_space.jl")
include("test_models.jl")
include("test_chebyshev.jl")
include("test_dmt.jl")
