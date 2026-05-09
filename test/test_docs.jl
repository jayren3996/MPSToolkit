function _user_facing_documentation_files()
  root = normpath(joinpath(@__DIR__, ".."))
  files = [
    joinpath(root, "docs", "Manifest.toml"),
    joinpath(root, "docs", "index.md"),
    joinpath(root, "docs", "api.md"),
    joinpath(root, "docs", "architecture.md"),
    joinpath(root, "docs", "examples.md"),
  ]

  for (dir, _, names) in walkdir(joinpath(root, "docs", "src"))
    for name in names
      endswith(name, ".md") && push!(files, joinpath(dir, name))
    end
  end

  return root, files
end

@testset "user-facing documentation is portable" begin
  root, files = _user_facing_documentation_files()

  for file in files
    text = read(file, String)

    @testset "$(relpath(file, root))" begin
      @test findfirst("/Users/", text) === nothing
    end
  end
end
