using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random
using Test

@testset "DMT checkpoint/restart" begin
  # A run at chi >= 400, d >= 4 takes hours, and a cluster wall-clock limit that kills it at
  # hour eight with nothing on disk is the failure mode these exist for. The interesting part
  # is not the I/O -- it is refusing to resume across a parameter change.
  z3 = Matrix{ComplexF64}(Diagonal([1.0, 0.0, -1.0]))
  function sample_state(nsites = 6)
    Random.seed!(4242)
    sites = operator_siteinds(nsites; d = 3)
    return add(operator_basis_state(sites, fill(1, nsites)),
      0.3 * random_mps(sites; linkdims = 8); maxdim = 8, cutoff = 0.0)
  end
  parameters = Dict(:nsites => 6, :d => 3, :maxdim => 40, :dt => 0.05)

  @testset "round trip preserves the state exactly" begin
    dir = mktempdir()
    rho = sample_state()
    before = operator_expectation_profile(rho, [(2, z3)]; normalize = false)[1]
    path = dmt_checkpoint_save(dmt_checkpoint_path("run", 1.2; dir = dir), rho, parameters;
      time = 1.2, sweep = 6, observables = Dict(:trace => [1.0, 1.0]))
    @test isfile(path)
    restored = dmt_checkpoint_load(path)
    @test restored.time == 1.2
    @test restored.sweep == 6
    @test restored.observables[:trace] == [1.0, 1.0]
    @test maxlinkdim(restored.state) == maxlinkdim(rho)
    # The element type has to survive too, or a resumed run silently drops off the real
    # arithmetic path and runs ~2x slower than the one that wrote the checkpoint.
    @test mapreduce(eltype, promote_type, ITensorMPS.data(restored.state)) === Float64
    # The restored state must be usable, not merely equal: its site indices live inside its own
    # tensors, and a fresh operator_siteinds call would mint ids that will not contract.
    after = operator_expectation_profile(restored.state, [(2, z3)]; normalize = false)[1]
    @test after ≈ before rtol = 1e-14
  end

  @testset "run cadence and provenance survive a restart" begin
    dir = mktempdir()
    rho = sample_state()
    run_state = ConstrainedDMTRunState(
      completed_steps=7, physical_time=0.7, steps_since_projection=1,
      projection_count=2)
    provenance = dmt_run_provenance(algorithm_version="pxp-dmt-v1")
    path = dmt_checkpoint_save(dmt_checkpoint_path("stateful", 0.7; dir=dir), rho,
      parameters; time=0.7, sweep=7, run_state=run_state, provenance=provenance)
    restored = dmt_checkpoint_load(path)
    restored_state = dmt_checkpoint_run_state(restored)
    @test restored_state.completed_steps == 7
    @test restored_state.physical_time == 0.7
    @test restored_state.steps_since_projection == 1
    @test restored_state.projection_count == 2
    @test dmt_checkpoint_provenance(restored)[:algorithm_version] == "pxp-dmt-v1"
    @test dmt_checkpoint_provenance(restored)[:code_sha] isa Union{Nothing,String}
    @test_throws ArgumentError dmt_checkpoint_save(
      dmt_checkpoint_path("bad-state", 0.8; dir=dir), rho, parameters;
      time=0.8, sweep=8, run_state=run_state)

    # Public data can be carried through a complete load -> advance -> save -> load cycle
    # without exposing or accidentally resubmitting the reserved metadata keys.
    carried = dmt_checkpoint_observables(restored)
    carried[:trace] = [1.0]
    restored_state.completed_steps += 1
    restored_state.physical_time += 0.1
    second_path = dmt_checkpoint_save(dmt_checkpoint_path("stateful", 0.8; dir=dir),
      restored.state, parameters; time=0.8, sweep=8, observables=carried,
      run_state=restored_state, provenance=dmt_checkpoint_provenance(restored))
    second = dmt_checkpoint_load(second_path)
    @test dmt_checkpoint_observables(second) == Dict(:trace => [1.0])
    @test dmt_checkpoint_run_state(second).completed_steps == 8
    @test dmt_checkpoint_provenance(second)[:algorithm_version] == "pxp-dmt-v1"

    legacy_path = dmt_checkpoint_save(dmt_checkpoint_path("legacy", 0.0; dir=dir), rho,
      parameters; time=0.0, sweep=0)
    legacy = dmt_checkpoint_load(legacy_path)
    @test dmt_checkpoint_run_state(legacy) === nothing
    @test isempty(dmt_checkpoint_provenance(legacy))
  end

  @testset "resume finds the newest checkpoint" begin
    dir = mktempdir()
    rho = sample_state()
    @test dmt_checkpoint_resume("run", parameters; dir = dir) === nothing
    for (time, sweep) in ((0.4, 2), (2.4, 12), (1.2, 6))     # deliberately out of order
      dmt_checkpoint_save(dmt_checkpoint_path("run", time; dir = dir), rho, parameters;
        time = time, sweep = sweep)
    end
    # Times are written fixed-width so plain string ordering is chronological; a naive format
    # would make "t_10.0" sort before "t_2.4" and silently resume from the wrong snapshot.
    @test dmt_checkpoint_resume("run", parameters; dir = dir).time == 2.4
    @test dmt_checkpoint_resume("other_run", parameters; dir = dir) === nothing
  end

  @testset "a parameter change is refused, not silently resumed" begin
    dir = mktempdir()
    dmt_checkpoint_save(dmt_checkpoint_path("run", 1.2; dir = dir), sample_state(), parameters;
      time = 1.2, sweep = 6)
    for changed in (Dict(:nsites => 6, :d => 3, :maxdim => 80, :dt => 0.05),   # different budget
      Dict(:nsites => 8, :d => 3, :maxdim => 40, :dt => 0.05),                 # different chain
      Dict(:nsites => 6, :d => 3, :maxdim => 40))                              # dropped a key
      failure = try
        dmt_checkpoint_resume("run", changed; dir = dir)
        nothing
      catch exception
        exception
      end
      @test failure isa ArgumentError
      # The message must name what differs: the whole point is telling someone which knob they
      # moved, since the alternative is a time series with a discontinuity nothing can detect.
      @test occursin("written by a different run", sprint(showerror, failure))
    end
  end

  @testset "a half-written checkpoint is never mistaken for a good one" begin
    dir = mktempdir()
    dmt_checkpoint_save(dmt_checkpoint_path("run", 1.2; dir = dir), sample_state(), parameters;
      time = 1.2, sweep = 6)
    # The write goes to a scratch name and is moved into place, so a kill mid-write leaves a
    # .partial rather than a truncated .dmt. Nothing should be left behind on success.
    @test isempty(filter(f -> endswith(f, ".partial"), readdir(dir)))
    @test length(filter(f -> endswith(f, ".dmt"), readdir(dir))) == 1
    # A file that is not a checkpoint at all is rejected rather than deserialized blindly.
    junk = joinpath(dir, "junk_t_00000.000.dmt")
    write(junk, "not a checkpoint")
    @test_throws Exception dmt_checkpoint_load(junk)
  end
end
